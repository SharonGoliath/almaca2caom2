# -*- coding: utf-8 -*-
# ***********************************************************************
# ******************  CANADIAN ASTRONOMY DATA CENTRE  *******************
# *************  CENTRE CANADIEN DE DONNÉES ASTRONOMIQUES  **************
#
#  (c) 2018.                            (c) 2018.
#  Government of Canada                 Gouvernement du Canada
#  National Research Council            Conseil national de recherches
#  Ottawa, Canada, K1A 0R6              Ottawa, Canada, K1A 0R6
#  All rights reserved                  Tous droits réservés
#
#  NRC disclaims any warranties,        Le CNRC dénie toute garantie
#  expressed, implied, or               énoncée, implicite ou légale,
#  statutory, of any kind with          de quelque nature que ce
#  respect to the software,             soit, concernant le logiciel,
#  including without limitation         y compris sans restriction
#  any warranty of merchantability      toute garantie de valeur
#  or fitness for a particular          marchande ou de pertinence
#  purpose. NRC shall not be            pour un usage particulier.
#  liable in any event for any          Le CNRC ne pourra en aucun cas
#  damages, whether direct or           être tenu responsable de tout
#  indirect, special or general,        dommage, direct ou indirect,
#  consequential or incidental,         particulier ou général,
#  arising from the use of the          accessoire ou fortuit, résultant
#  software.  Neither the name          de l'utilisation du logiciel. Ni
#  of the National Research             le nom du Conseil National de
#  Council of Canada nor the            Recherches du Canada ni les noms
#  names of its contributors may        de ses  participants ne peuvent
#  be used to endorse or promote        être utilisés pour approuver ou
#  products derived from this           promouvoir les produits dérivés
#  software without specific prior      de ce logiciel sans autorisation
#  written permission.                  préalable et particulière
#                                       par écrit.
#
#  This file is part of the             Ce fichier fait partie du projet
#  OpenCADC project.                    OpenCADC.
#
#  OpenCADC is free software:           OpenCADC est un logiciel libre ;
#  you can redistribute it and/or       vous pouvez le redistribuer ou le
#  modify it under the terms of         modifier suivant les termes de
#  the GNU Affero General Public        la “GNU Affero General Public
#  License as published by the          License” telle que publiée
#  Free Software Foundation,            par la Free Software Foundation
#  either version 3 of the              : soit la version 3 de cette
#  License, or (at your option)         licence, soit (à votre gré)
#  any later version.                   toute version ultérieure.
#
#  OpenCADC is distributed in the       OpenCADC est distribué
#  hope that it will be useful,         dans l’espoir qu’il vous
#  but WITHOUT ANY WARRANTY;            sera utile, mais SANS AUCUNE
#  without even the implied             GARANTIE : sans même la garantie
#  warranty of MERCHANTABILITY          implicite de COMMERCIALISABILITÉ
#  or FITNESS FOR A PARTICULAR          ni d’ADÉQUATION À UN OBJECTIF
#  PURPOSE.  See the GNU Affero         PARTICULIER. Consultez la Licence
#  General Public License for           Générale Publique GNU Affero
#  more details.                        pour plus de détails.
#
#  You should have received             Vous devriez avoir reçu une
#  a copy of the GNU Affero             copie de la Licence Générale
#  General Public License along         Publique GNU Affero avec
#  with OpenCADC.  If not, see          OpenCADC ; si ce n’est
#  <http://www.gnu.org/licenses/>.      pas le cas, consultez :
#                                       <http://www.gnu.org/licenses/>.
#
#  $Revision: 4 $
#
# ***********************************************************************
#

import logging
import numpy
import os

from datetime import datetime

from astropy import units, time, constants
from caom2 import Energy, EnergyBand, shape, Interval, Algorithm
from caom2 import Telescope, Instrument, Target, Proposal, Position
from caom2 import DerivedObservation, Provenance, Artifact, Plane
from caom2 import DataProductType, CalibrationLevel, TargetType
from caom2 import ReleaseType, ProductType, ObservationIntentType
from caom2 import Time, ObservationURI, Environment, PlaneURI

from caom2pipe import astro_composable as ac
from caom2pipe import manage_composable as mc

ARCHIVE = 'ALMACA'


class AlmacaName(mc.StorageName):
    """
    Naming rules:
    - product id should be the file name without the “.ms.split.cal” suffix,
    and the uid___ prefix.

    ALMA IDs:
    - science project id
    - science goal ous id
    - group ous id - GOUS - contains all MOUS needed to achieve a science goal
        - can include different configurations or arrays
    - member ous id - MPUS - set of all scheduling blocks needed to achieve
        part of a science goal, usually associated with a specific array
        configuration
    - ASDM id - ALMA Science Data Model - final product of each observation
        in the archive
    - OUS Obs Unit Set

    - if calibrate science observation of a MOUS have to be combined with
      science observations of another MOUS, they are grouped together in a
      GOUS

    # To harmonize with 'ALMA' collection:
    #     - observation_id value is the ALMA Member ous id.
    #
    #     - product_name value is the ALMA Member ous id, plus the spot
    #       in the MS name occupied by values like 'highres_spw2'
    #
    #     -

    HK - 04-02-20
    #
    # On Observation cardinality:

    A002_Xb999fd_X602.CAL.J1751+0939.highres_spw2 and
    A002_Xb999fd_X602.CAL.J1751+0939.highres_spw3 are part of the same
    observation, but neither A002_Xb999fd_X602.CAL.J1924-2914.highres_spw2 or
    A002_Xb999fd_X602.SCI.J1851+0035.highres_spw2 would be.
    """
    def __init__(self, fname_on_disk=None, file_name=None, obs_id=None,
                 file_id=None):
        super(AlmacaName, self).__init__(collection=ARCHIVE,
                                         collection_pattern='*',
                                         fname_on_disk=fname_on_disk)
        self._file_name = os.path.basename(fname_on_disk)
        temp = self._file_name.split('.')
        if len(temp) < 5:
            raise mc.CadcException('Not a split product.')
        asdm_str = temp[0].replace('uid___', '')
        self._obs_id = f'{asdm_str}.{temp[1]}.{temp[2]}'

        # TODO - hard-coded for single-band splitting testing right now
        self._science_goal_id = 'uid://A001/X88b/X21'
        self._group_id = 'uid://A001/X88b/X22'
        self._mous_id = 'uid://A001/X88b/X23'
        self._asdm_id = f'uid://{asdm_str.replace("_", "/")}'

        self._product_id = temp[3]
        self._intent = (ProductType.CALIBRATION if '.CAL.' in fname_on_disk
                        else ProductType.SCIENCE)
        self._ms = fname_on_disk
        self._log_dir = '/data/calibrated'
        self._input_ms_metadata = f'{self._log_dir}/{temp[0]}' \
                                  f'.ms.split.cal/md.pk'
        self._logger = logging.getLogger(__name__)
        self._logger.error(self)

    def __str__(self):
        return f'\nobs_id {self.obs_id}' \
               f'\nfile_name {self.file_name}'

    @property
    def file_name(self):
        return self._file_name

    @property
    def input_ms_metadata(self):
        return self._input_ms_metadata

    @property
    def input_uri(self):
        return PlaneURI(mc.CaomName.make_plane_uri(
            ARCHIVE, self._obs_id, self._product_id))

    @property
    def intent(self):
        return self._intent

    @property
    def log_dir(self):
        return self._log_dir

    @property
    def mous_id(self):
        return self._mous_id

    @property
    def product_id(self):
        return self._product_id

    @property
    def uri(self):
        return mc.build_uri(ARCHIVE, '{}.tar.gz'.format(self._ms))

    @staticmethod
    def to_uid(value):
        return value.replace('___', '://').replace('_', '/', 2)


def _get_band_name(override):
    return 'Band {}'.format(override.get('band'))


def build_position(db_content, field_index, md_name):
    f_names = os.listdir(os.path.dirname(md_name))
    index = -1
    for f_name in f_names:
        if 'listobs' in f_name:
            fqn = '{}/{}'.format(os.path.dirname(md_name), f_name)
            with open(fqn) as f:
                temp = f.readlines()
            for ii in temp:
                if 'RA               Decl' in ii:
                    index = temp.index(ii) + 1
                    break
            break

    if index == -1:
        raise mc.CadcException('ra dec not found')

    ra = temp[index].split()[3]
    dec = temp[index].split()[4].replace('.', ':', 2)
    result_ra, result_dec = ac.build_ra_dec_as_deg(ra, dec)

    # HK 19-08-19
    # Looking at the ALMA web archive listing, it claims a 'FOV' (field of
    # view) of 53.38 arcsec, which would presumably be the diameter of the
    # observed area. That corresponds to a radius of 26.69, which is a
    # little larger than the value of 24 given in caom2.

    fov = db_content['Field of view'][field_index]
    radius = ((fov / 2.0) * units.arcsec).to(units.degree)
    bounds = shape.Circle(
        center=shape.Point(result_ra, result_dec),
        radius=radius.value)

    # HK 17-10-19
    # Position: resolution: should not be null.  For the sample dataset that
    # we've been iterating on, the alma web query lists a value of 0.4
    # arcsec (the 'Ang. res.' parameter).  While this number in principle
    # varies a little bit for the different spectral windows, I am not seeing
    # an easy way to extract more precise values using msmd.  Since the
    # achievable angular resolution actually depends a bit on how different
    # baselines are weighted during imaging, there isn't an exact fixed value
    # no matter what, so I think it would be fine to just use the ALMA web
    # query value as being 'good enough'.
    resolution = db_content['Spatial resolution'][field_index]

    return Position(bounds=bounds,
                    sample_size=None,
                    time_dependent=False,
                    resolution=resolution)


def build_energy(override):

    spectral_windows = override.get('spectral_windows')
    sample_size = mc.to_float(override.get('energy_sample_size'))

    # HK 19-08-19
    # I'm still not quite sure that I follow this one.  I understand your
    # point from earlier about merging together overlapping wavelength
    # ranges, and that should be fine, although it might hide some
    # important information in the energy:resolution parameter, since
    # each of the overlapping ranges may have different spectral
    # resolutions.  By my approximate calculations, the 4 wavelength
    # ranges covered are 0.00259 to 0.00263, 0.00263 to 0.00267, 0.0026025
    # to 0.0026038, and 0.0026018 to 0.0026044.  If I were to merge those,
    # I would get 0.00259 to 0.00267.  I don't understand how the caom2
    # model lists 0.00259 to 0.002602 and then 0.002626 and 0.002672.
    # Specifically, the wavelengths covered by the first of the 4 ranges
    # I list, already span a much larger range than the first range listed
    # in caom2.

    energy = Energy()
    energy.em_band = EnergyBand.MILLIMETER
    energy.dimension = 1

    wvlns = []
    mid_wvln = []

    for spw in spectral_windows:
        wvln = numpy.array((_from_hz_to_m(spw[0]), _from_hz_to_m(spw[1])))
        wvlns.append(wvln)
        mid_wvln.append(wvln[0] + wvln[1])
    order = numpy.argsort(mid_wvln)

    min_bound = None
    max_bound = None
    si = []
    for idx in order:
        lower = min(wvlns[idx])
        upper = max(wvlns[idx])

        si = _add_subinterval(si, (lower, upper))
        if min_bound is not None:
            min_bound = min(min_bound, lower)
        else:
            min_bound = lower
        if max_bound is not None:
            max_bound = max(max_bound, upper)
        else:
            max_bound = upper

    samples = []
    for s in si:
        samples.append(shape.SubInterval(s[0], s[1]))

    energy.bounds = Interval(min_bound, max_bound, samples=samples)

    # HK 15-10-19
    #
    # It looks like the caom2 model puts energy in wavelength units (I'm not
    # clear whether that is metres or centimetres?), whereas the numbers
    # I quoted above are directly from msmd and are given in frequency units
    # of Hertz.  I'm not sure what level of precision you would use for the
    # conversion, but here's roughly what you'd want to do:
    #
    # [resolution in wavelength units] / [central wavelength] =
    # [resolution in frequency units] /[central frequency]
    #
    # Using [central wavelength] = [speed of light] / [central frequency], and
    # a speed of light of 2.9979e8 m/s (or whatever precision you need, and
    # convert to cm/s if needed), you should be able to run the calculation
    # with values you've already extracted.  For a central frequency of
    # 113GHz, I get a value of about 3.7e-7 m, assuming I've done my quick
    # calculation correctly.
    #
    # NB: since both the sampleSize and resolution parameters are looking at a
    # differential wavelength measurement, both conversions would follow the
    # same formula as I've written.
    #
    mean_frequency = (_from_m_to_hz(min_bound) + _from_m_to_hz(max_bound)) / 2
    energy.sample_size = _delta_hz_to_m(sample_size, mean_frequency)
    energy_resolution = mc.to_float(override.get('energy_resolution'))

    # HK 03-12-19
    # resolving power is unit-less
    # ResolvingPower = mean[frequency_Hz] / chanres_Hz
    energy.resolving_power = mean_frequency / energy_resolution

    # HK 3-10-19
    # energy: bandpassName: could this also be Band3?  (I know it is already
    # listed under 'instrument' in the top level plane)
    energy.bandpass_name = _get_band_name(override)
    return energy


def _delta_hz_to_m(from_value, mean_frequency):
    # HK - 09-12-19
    # Use the following:
    #
    # (change_in_frequency)/(frequency) = (change_in_wavelength)/(wavelength)
    #
    # You can convert that to the following:
    #
    # (change_in_wavelength) =
    #   (change_in_frequency)*(speed_of_light)/(frequency^2)
    #
    # change in frequency => 'energy_sample_size' => msmd.chanwidth
    # frequency => mean of 'spectral_windows' => msmd.chanfreqs
    return ((from_value * constants.c) / (mean_frequency ** 2)).value


def _from_hz_to_m(val):
    # JJK 15-10-19
    # I'll just interject into here that what the pyCAOM code should do (and
    # I really mean should) is accept values that have units when
    # instantiating the object and then convert that to the correct unit when
    # building the XML.   (turn away Helen!)  This could be achieved VERY
    # easily in python using astropy.units and would save a large number of
    # foot-gun situations.  Eg.
    #
    # caom2.plane.Energy(bounds=[lower_frequency*units.GHz,
    #                    upper_frequency*units.GHz])
    #
    # and then Energy uses those values as lower_frequency.to(units.m).value
    # when sending to the service.
    #
    # If a value (rather than a Quantity) is sent to Energy (ie no units)
    # then Energy can assume they are in m or through an error, but for
    # backwards friendly it should accept this as unit.m by default and cast
    # the value to a unit.m Quantity
    return (val * units.Hz).to(units.m, equivalencies=units.spectral()).value


def _from_m_to_hz(val):
    return (val * units.m).to(units.Hz, equivalencies=units.spectral()).value


def _add_subinterval(si_list, subinterval):
    # Adds an interval to a list of intervals eliminating (merging) any
    # overlaps

    if not si_list:
        return [subinterval]
    # check for overlaps
    # beginning of the list?
    if subinterval[1] < si_list[0][0]:
        return [subinterval] + si_list
    if subinterval[0] > si_list[-1][1]:
        return si_list + [subinterval]
    result = []
    for si in si_list:
        if (subinterval[0] <= si[0] <= subinterval[1]) or \
                (si[0] <= subinterval[0] <= si[1]):
            # overlap detected
            subinterval = (min(si[0], subinterval[0]),
                           max(si[1], subinterval[1]))
        else:
            if subinterval[0] < si[0]:
                result += [subinterval]
                result += si_list[si_list.index(si):]
                return result
            else:
                result += [si]
    return result + [subinterval]


def read_md_pk(fqn):
    temp = mc.read_from_file(fqn)
    assert temp is not None, 'expected result'
    result = {'spectral_windows': []}
    for line in temp:
        temp2 = line.split(',', 1)
        if temp2[0].strip() == 'spectral_windows':
            x = temp2[1].strip().split(',')
            count = 0
            while count < len(x):
                y = (mc.to_float(x[count].replace('(', '').replace('[', '').replace(']', '').replace(')', '')),
                     mc.to_float(x[count+1].replace('(', '').replace('[', '').replace(']', '').replace(')', '')))
                count += 2
                result[temp2[0]].append(y)
        elif type(temp2[1]) is str:
            result[temp2[0]] = temp2[1].strip()
        else:
            result[temp2[0]] = temp2[1]
    return result


def build_time(override, almaca_name):

    # HK 05-09-19
    # For the time dimension:
    # - resolution: this is far less clear to me, but I think in principle,
    # once could try to make an image using some time-based subset of the full
    # measurement set.  So in principle, I suppose one could make [sampleSize]
    # independent images.  In practice, there probably would be insufficient
    # data to actually make decent images for such a small subset of the
    # data.  Something like the number of times that the source is observed in
    # between calibrators would probably be a more appropriate / practical
    # level of time sampling, but I don't think there would be an easy way to
    # pull that information out of the metadata.  The time resolution value
    # could be listed as 'null' to avoid the confusion.
    resolution = None

    # HK 04-02-20
    # Use time values from the input MS.
    input_meta_data = read_md_pk(almaca_name.input_ms_metadata)
    start_date = mc.to_float(input_meta_data.get('start_date'))
    end_date = mc.to_float(input_meta_data.get('end_date'))

    # HK 14-08-09
    # If 'exposure' is supposed to be the total (useful) exposure time
    # on-source, then this parameter should be msmd.effexposuretime().
    # The 'itime' parameter calculated in get_msmd.py gives a larger
    # value, as I believe it includes the full time on-source, regardless
    # of what mode the data is being taken in.
    exposure_time = mc.to_float(override.get('effexposuretime'))

    time_bounds = Interval(start_date,
                           end_date,
                           samples=[shape.SubInterval(start_date, end_date)])
    return Time(bounds=time_bounds,
                dimension=1,
                resolution=resolution,
                sample_size=mc.to_float(override.get('time_sample_size')),
                exposure=exposure_time)


def _build_obs(override, db_content, fqn, index, almaca_name):

    obs_date = db_content['Observation date'][index]
    if obs_date is None:
        raise mc.CadcException('No observation date for {}'.format(fqn))
    else:
        obs_date = time.Time(obs_date).to_datetime()

    # of_site('alma')
    # 2225015.30883296, -5440016.41799762, -2481631.27428014
    #
    size = db_content['Array'][index]
    telescope = Telescope(name="ALMA-{}".format(size),
                          geo_location_x=2225142.18,
                          geo_location_y=-5440307.37,
                          geo_location_z=-2481029.852)

    instrument = Instrument(name=_get_band_name(override))

    # HK - 14-08-19
    # target: this should be the science target / science fieldname and
    # not the calibrator fieldname, as it is presently (i.e., should be
    # J1851+0035, not J1924-2914).  Or, we may need to continue to leave
    # that field blank at this level.  A single calibrated measurement set
    # may contain multiple science target names, which would not be
    # properly captured at this level.  [For the raw data, the target
    # field is left blank]
    target = Target(name=override.get('field'),
                    standard=False,
                    moving=False,
                    target_type=TargetType.OBJECT)

    # db_content as votable:
    # >>> t.colnames
    # ['Project_code', 'Source_name', 'RA', 'Dec', 'Galactic_longitude',
    # 'Galactic_latitude', 'Band', 'Spatial_resolution',
    # 'Frequency_resolution', 'Array', 'Mosaic', 'Integration',
    # 'Release_date', 'Frequency_support', 'Velocity_resolution',
    # 'Pol_products', 'Observation_date', 'PI_name', 'SB_name',
    # 'Proposal_authors', 'Line_sensitivity__10_km_s_',
    # 'Continuum_sensitivity', 'PWV', 'Group_ous_id', 'Member_ous_id',
    # 'Asdm_uid', 'Project_title', 'Project_type', 'Scan_intent',
    # 'Field_of_view', 'Largest_angular_scale', 'QA2_Status', 'COUNT',
    # 'Science_keyword', 'Scientific_category', 'ASA_PROJECT_CODE']
    #
    # db_content.colnames as html:
    #
    # ['Project code', 'Source name', 'RA', 'Dec', 'Galactic longitude',
    # 'Galactic latitude', 'Band', 'Spatial resolution',
    # 'Frequency resolution', 'Array', 'Mosaic', 'Integration',
    # 'Release date', 'Frequency support', 'Velocity resolution',
    # 'Pol products', 'Observation date', 'PI name', 'SB name',
    # 'Proposal authors', 'Line sensitivity (10 km/s)',
    # 'Continuum sensitivity', 'PWV', 'Group ous id', 'Member ous id',
    # 'Asdm uid', 'Project title', 'Project type', 'Scan intent',
    # 'Field of view', 'Largest angular scale', 'QA2 Status', 'COUNT',
    # 'Science keyword', 'Scientific category', 'ASA_PROJECT_CODE']

    # HK - 14-08-19
    # can we include the project code, and not just the observation UID
    # somewhere in here?  For the raw data, it looks like the project
    # code was included as 'proposal: ID', whereas for the calibrated
    # measurement set, proposal: ID is now set to the UID and the project
    # code (2016.1.00010.S) is not captured anywhere.  Could the 'project'
    # field, currently set as 'null' be used for this?
    proposal = Proposal(id=db_content['Project code'][index],
                        project=override.get('project'),
                        pi_name=db_content['PI name'][index],
                        title=db_content['Project title'][index])

    keywords = db_content['Science keyword'][index]
    if keywords is not None:
        proposal.keywords = set(keywords.split())

    environment = Environment()
    environment.tau = db_content['PWV'][index]/0.935 + 0.35
    environment.wavelength_tau = 350*units.um.to(units.meter)

    intent = (ObservationIntentType.SCIENCE if almaca_name.intent is
              ProductType.SCIENCE else ObservationIntentType.CALIBRATION)

    algorithm = Algorithm(name='single band split')
    #
    # PD, SG 15-08-19
    # make it a composite, algorithm name something like
    # 'target splitting'
    #
    observation = DerivedObservation(collection=ARCHIVE,
                                     observation_id=almaca_name.obs_id,
                                     sequence_number=None,
                                     intent=intent,
                                     type="OBJECT",
                                     proposal=proposal,
                                     telescope=telescope,
                                     instrument=instrument,
                                     target=target,
                                     meta_release=obs_date,
                                     algorithm=algorithm,
                                     environment=environment)
    observation.members.add(
        ObservationURI(
            mc.CaomName.make_obs_uri_from_obs_id('ALMA', 'A001_X88b_X23')))
    return observation


def get_provenance(almaca_name):
    # HK 14-08-19
    # provenance: version - capture the information on what version of
    # CASA was used to run the calibration script.  We might appreciate
    # having that information saved later on (as might an advanced user).
    # This would be possible to capture from the 'casa[date].log' file
    # generated automatically during processing - the second line
    # includes 'CASA version XXX'.
    version_result = None
    last_result = None
    log_dir = almaca_name.log_dir
    logging.error('checking {}'.format(log_dir))
    if os.path.exists(log_dir):
        # logging.error('exists {}'.format(log_dir))
        log_dir_contents = os.listdir(log_dir)
        for ii in log_dir_contents:
            if ii.startswith('casa-') and ii.endswith('.log'):
                log_fqn = '{}/{}'.format(log_dir, ii)
                if os.path.exists(log_fqn):
                    with open(log_fqn, 'r') as f:
                        temp = f.readlines()
                    for jj in temp:
                        if 'CASA Version' in jj:
                            version_result = jj.split('CASA Version ')[1]

                    # get the timestamp from the filename, use it as the
                    # 'last_executed'
                    temp = ii.replace('casa-', '').replace('.log', '')
                    last_result = datetime.fromtimestamp(mc.make_seconds(temp))
    # TODO time.Time(override.get('casa_run_date')).datetime

    # The rest of the MAG seemed less concerned about the various OUS IDs being
    # searchable within the archive.  I think it would still be best to include
    # the information somewhere just in case.  My guess is that the ASDM UID is
    # the most important one to be searchable, and that it would also be quite
    # appropriate to be listed as the 'reference' under 'provenance'.  (It might
    # even eventually be linked directly to the associated raw data file.)  The
    # rest of the science/group/member OUS IDs could perhaps be listed within
    # the keywords section like this:
    #
    # ScienceGoalOUSID: [ugly string]; GroupOUSID: [ugly string#2];
    # MemberOUSID: [ugly string#3]  (or whatever formatting will work within
    # the keyword field).

    provenance = Provenance(name='CASA',
                            version=version_result,
                            last_executed=last_result,
                            reference='https://casa.nrao.edu/')
    provenance.keywords.add(f'ScienceGoalOUSID: {almaca_name._science_goal_id}')
    provenance.keywords.add(f'GroupOUSID: {almaca_name._group_id}')
    provenance.keywords.add(f'MemberOUSID: {almaca_name._mous_id}')
    provenance.keywords.add(f'ASDM ID: {almaca_name._asdm_id}')
    return provenance


def _get_index(almaca_name, db_content):
    # logging.error(f'almaca_name is {almaca_name}')
    # almaca_name is obs_id A002_Xb999fd_X602.CAL.J1924-2914,
    # fname_on_disk /data/uid___A002_Xb999fd_X602.CAL.J1924-2914.highres_spw2.
    # ms.split.cal,
    # file_name A002_Xb999fd_X602.CAL.J1924-2914.fits,
    # lineage highres_spw2/ad:ALMACA/A002_Xb999fd_X602.CAL.J1924-2914.fits.gz

    # logging.error(f'almaca_name.mous_id {almaca_name.mous_id}')
    # uid://A002/Xb999fd/X602.CAL.J1924-2914.highres/spw2.ms.split.cal

    # logging.error(db_content['Member ous id'])
    # -------------------
    # uid://A001/X88b/X37
    # uid://A001/X88b/X2b
    # uid://A001/X88b/X27
    # uid://A001/X88b/X33
    # uid://A001/X88b/X23
    # uid://A001/X88b/X2f

    # HK - conversation - 10-09-19
    # use the Member_ous_id for the field index
    index = None
    count = 0
    found = False
    for ii in db_content['Member ous id']:
        if ii == almaca_name.mous_id:
            index = count
            found = True
            break
        count += 1

    if not found:
        raise mc.CadcException(
            'Could not find field index for {}'.format(almaca_name.mous_id))
    return index


def build_observation(db_content, observation, md_name):

    override = read_md_pk(md_name)

    fqn = override.get('fqn')
    almaca_name = AlmacaName(fname_on_disk=fqn)
    # logging.error(db_content.colnames)
    # logging.error('fqn is {}'.format(fqn))
    field_index = _get_index(almaca_name, db_content)
    # field_index = 0
    if observation is None:
        observation = _build_obs(override, db_content, fqn, field_index,
                                 almaca_name)

    provenance = get_provenance(almaca_name)
    provenance.inputs.add(PlaneURI('caom:ALMA/A001_X88b_X23/A001_X88b_X23-raw'))

    release_date = db_content['Release date'][field_index]
    if release_date is None:
        raise mc.CadcException('No release date for {}'.format(fqn))
    else:
        release_date = time.Time(release_date).to_datetime()

    logging.error('Add plane {} to {}'.format(almaca_name.product_id,
                                              almaca_name.obs_id))
    plane = Plane(product_id=almaca_name.product_id,
                  data_release=release_date,
                  meta_release=observation.meta_release,
                  provenance=provenance)

    plane.position = build_position(db_content, field_index, md_name)
    plane.energy = build_energy(override)
    plane.polarization = None
    plane.time = build_time(override, almaca_name)

    # HK 14-08-2019
    # dataProductType should be 'visibility'
    plane.data_product_type = DataProductType.VISIBILITY
    plane.calibration_level = CalibrationLevel.CALIBRATED

    observation.planes.add(plane)
    # TODO hard-coded
    observation.members.add(ObservationURI('caom:ALMA/A001_X88b_X23'))

    # HK 29-07-19
    # qa/ contains images, plots, and web page status views generated
    # during the original (non-CANFAR) calibration of the raw data.
    # We may want to consider retaining these files as well, as they give
    # a more advanced user an easier way to check on data quality,
    # potential issues with calibration, etc.  I believe they come
    # packaged with the rest of the 'products' tarball on the archive,
    # so they would be obtainable even if we do not keep a copy.  These
    # files are fairly small.
    # TODO override.get('artifact_uri')
    artifact = Artifact(
        uri=almaca_name.uri,
        product_type=almaca_name.intent,
        release_type=ReleaseType.DATA,
        content_type='application/x-tar',
        content_length=None)
    plane.artifacts.add(artifact)
    return observation
