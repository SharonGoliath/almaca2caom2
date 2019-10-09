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
import math
import numpy
import os

from datetime import datetime

from astropy import units, time, constants
from caom2 import SegmentType, Vertex, Point, Position, Polygon, shape
from caom2 import Energy, EnergyBand, shape, Interval, Algorithm
from caom2 import Telescope, Instrument, Target, Proposal
from caom2 import CompositeObservation, Provenance, Artifact, Plane
from caom2 import DataProductType, CalibrationLevel, TargetType
from caom2 import ReleaseType, ProductType, ObservationIntentType
from caom2 import Time, ObservationURI, Environment, PlaneURI

from caom2pipe import execute_composable as ec
from caom2pipe import manage_composable as mc

ARCHIVE = 'ALMACA'


class AlmacaName(ec.StorageName):
    """
    Naming rules:
    - product id should be the file name without the “.ms.split.cal” suffix,
    and the uid___ prefix.

    """
    def __init__(self, fname_on_disk=None, file_name=None, obs_id=None,
                 file_id=None):
        super(AlmacaName, self).__init__(collection=ARCHIVE,
                                         collection_pattern='*',
                                         fname_on_disk=fname_on_disk)
        temp = fname_on_disk.split('/')[-1].split('.')
        self._obs_id = '{}.{}.{}'.format(
            temp[0], temp[1], temp[2]).replace('uid___', '')
        self._product_id = temp[3]
        self._mous_id = AlmacaName.to_uid(
            fname_on_disk.split('/calibrated')[0].split('/')[-1].replace(
                'member.', ''))

    @property
    def product_id(self):
        return self._product_id

    @property
    def input_uri(self):
        return PlaneURI(ec.CaomName.make_plane_uri(
            ARCHIVE, self._obs_id, self._product_id))

    @property
    def mous_id(self):
        return self._mous_id

    @staticmethod
    def to_uid(value):
        return value.replace('___', '://').replace('_', '/')


def _get_band_name(override):
    return 'Band {}'.format(override.get('band'))


def build_position(fqn, db_content, field_index, md_name):
    # f_name = md_name.replace('md', 'listobs_').replace('.pk', '.txt')
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
    # logging.error('ra {} dec {}'.format(ra, dec))

    from astropy import units
    from astropy.coordinates import SkyCoord

    result = SkyCoord(ra, dec, frame='icrs',
                      unit=(units.hourangle, units.deg))

    # HK 19-08-19
    # Looking at the ALMA web archive listing, it claims a 'FOV' (field of
    # view) of 53.38 arcsec, which would presumably be the diameter of the
    # observed area. That corresponds to a radius of 26.69, which is a
    # little larger than the value of 24 given in caom2.

    # , radius=24 * units.arcsec
    fov = db_content['Field of view'][field_index]
    radius = (fov / 2.0) * units.arcsec
    bounds = shape.Circle(
        center=shape.Point(result.ra.degree, result.dec.degree),
        radius=radius.value)
    return Position(bounds=bounds,
                    sample_size=None,
                    time_dependent=False)


def build_energy(override):

    spectral_windows = override.get('spectral_windows')
    sample_size = override.get('energy_sample_size')

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
    c = constants.c.to('m/s').value

    wvlns = []
    mid_wvln = []

    for spw in spectral_windows:
        wvln = numpy.array((c / spw[0], c / spw[1]))
        logging.error('wvln {} spw {}'.format(wvln, spw))
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
    energy.sample_size = mc.to_float(sample_size)
    energy.resolving_power = mc.to_float(override.get('energy_resolution'))

    # HK 3-10-19
    # energy: bandpassName: could this also be Band3?  (I know it is already
    # listed under 'instrument' in the top level plane)
    energy.bandpass_name = _get_band_name(override)
    return energy


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


def build_time(override):

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

    start_date = mc.to_float(override.get('start_date'))
    end_date = mc.to_float(override.get('end_date'))

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
                sample_size=mc.to_float(override.get('sample_size')),
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
    environment.tau = db_content['PWV'][0]/0.935 + 0.35
    environment.wavelength_tau = 350*units.um.to(units.meter)

    fqn = override.get('fqn')
    if fqn is None or '.SCI.' in fqn:
        intent = ObservationIntentType.SCIENCE
    else:
        intent = ObservationIntentType.CALIBRATION

    algorithm = Algorithm(name='band splitting')
    #
    # PD, SG 15-08-19
    # make it a composite, algorithm name something like
    # 'target splitting'
    #
    observation = CompositeObservation(collection=ARCHIVE,
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
            ec.CaomName.make_obs_uri_from_obs_id('ALMA', 'A001_X88b_X23')))
    return observation


def get_provenance(fqn):
    # HK 14-08-19
    # provenance: version - capture the information on what version of
    # CASA was used to run the calibration script.  We might appreciate
    # having that information saved later on (as might an advanced user).
    # This would be possible to capture from the 'casa[date].log' file
    # generated automatically during processing - the second line
    # includes 'CASA version XXX'.
    version_result = None
    last_result = None
    log_dir = '{}script'.format(fqn.split('calibrated')[0])
    # logging.error('checking {}'.format(log_dir))
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
                            version_result = jj.split(
                                'casa')[1].replace(':', '').strip()

                    # get the timestamp from the filename, use it as the
                    # 'last_executed'
                    temp = ii.replace('casa-', '').replace('.log', '')
                    last_result = datetime.fromtimestamp(mc.make_seconds(temp))
    # TODO time.Time(override.get('casa_run_date')).datetime
    return Provenance(name='CASA',
                      version=version_result,
                      last_executed=last_result,
                      reference='https://casa.nrao.edu/')


def _get_index(almaca_name, db_content):
    # HK - conversation - 10-09-19
    # use the Member_ous_id for the field index
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


def build_observation(override, db_content, observation, md_name):

    fqn = override.get('fqn')
    almaca_name = AlmacaName(fname_on_disk=fqn)
    # logging.error(db_content.colnames)
    # logging.error('fqn is {}'.format(fqn))
    field_index = _get_index(almaca_name, db_content)
    if observation is None:
        observation = _build_obs(override, db_content, fqn, field_index,
                                 almaca_name)

    provenance = get_provenance(fqn)

    if fqn is None or '.SCI.' in fqn:
        product_type = ProductType.SCIENCE
        # ip = AlmacaName(override.get('provenance'))
        # provenance.inputs.add(ip.input_uri)
    else:
        product_type = ProductType.CALIBRATION

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

    plane.position = build_position(fqn, db_content, field_index, md_name)
    plane.energy = build_energy(override)
    plane.polarization = None
    plane.time = build_time(override)

    # HK 14-08-2019
    # dataProductType should be 'visibility'
    plane.data_product_type = DataProductType.VISIBILITY
    plane.calibration_level = CalibrationLevel.CALIBRATED

    observation.planes.add(plane)

    uri_base = '{}.tar.gz'.format(fqn.split('/')[-1])

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
        uri=mc.build_uri(ARCHIVE, uri_base),
        product_type=product_type,
        release_type=ReleaseType.DATA,
        content_type='application/tar+gzip',
        content_length=None)
    plane.artifacts.add(artifact)
    return observation


# t = Table.read('./alma_query.xml', format='votable')
# t.info
# t.colnames - case sensitive 'Asdm_uid'
# for ii in t.colnames:
#     print(ii)
#     print(t[ii])
