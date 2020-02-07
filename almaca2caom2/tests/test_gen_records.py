import glob
import os
import math

ALMA_QUERY_DIR = '/usr/src/app/almaca2caom2/almaca2caom2/tests/data'
TEST_DATA_DIR = '/data'


def run_gen():
    import logging
    for root, dirs, files in os.walk(TEST_DATA_DIR):
        # logging.error('root is {}'.format(root))
        for dir_name in dirs:
            # logging.error('dir_name is {}'.format(dir_name))
            if dir_name.endswith('ms.split.cal'):
                fqn = '{}/{}'.format(root, dir_name)
                pk_file = '{}/md.pk'.format(fqn)
                if os.path.exists(pk_file):
                    os.unlink(pk_file)
                logging.error('fqn is {}'.format(fqn))
                provenance = None
                if '.SCI.' in fqn:
                    temp = fqn.replace('.SCI.', '.CAL.')
                    if os.path.exists(temp):
                        provenance = temp
                try:
                    logging.error('Processing {}'.format(fqn))
                    result = get_info(fqn, provenance)
                    logging.error(result)
                    try_vishead(fqn)
                    assert result is not None, 'expected result'
                    # pk_file = '{}/md.pk'.format(TEST_DATA_DIR)
                    logging.error('writing {}'.format(pk_file))
                    with open(pk_file, 'w') as pkl:
                        # pickle.dump(result, pkl)
                        # pkl.writelines(result)
                        for key in result.keys():
                            pkl.write('{}, {}\n'.format(key, result[key]))
                        pkl.write('{}, {}'.format('fqn', fqn))
                except Exception as e:
                    logging.error('failed {}'.format(fqn))
                    import traceback
                    logging.error(traceback.format_exc(e))


def run_gen_specific():
    fqn = '/data/for_CADC/2016.1.00010.S/science_goal.uid___A001_X88b_X21/' \
          'group.uid___A001_X88b_X22/member.uid___A001_X88b_X23/calibrated/' \
          'uid___A002_Xb945f7_X1551.SCI.J1851+0035.line_spw3.ms.split.cal'
    import logging
    pk_file = '{}/md.pk'.format(fqn)
    if os.path.exists(pk_file):
        os.unlink(pk_file)
    logging.error('fqn is {}'.format(fqn))
    provenance = None
    if '.SCI.' in fqn:
        temp = fqn.replace('.SCI.', '.CAL.')
        if os.path.exists(temp):
            provenance = temp
    try:
        logging.error('Processing {}'.format(fqn))
        result = get_info(fqn, provenance)
        logging.error(result)
        try_vishead(fqn)
        assert result is not None, 'expected result'
        # pk_file = '{}/md.pk'.format(TEST_DATA_DIR)
        logging.error('writing {}'.format(pk_file))
        with open(pk_file, 'w') as pkl:
            # pickle.dump(result, pkl)
            # pkl.writelines(result)
            for key in result.keys():
                pkl.write('{}, {}\n'.format(key, result[key]))
            pkl.write('{}, {}'.format('fqn', fqn))
    except Exception as e:
        logging.error('failed {}'.format(fqn))
        import traceback
        logging.error(traceback.format_exc(e))


def try_vishead(fqn):
    from datetime import datetime
    list_fqn = '{}/listobs.txt'.format(fqn)
    if os.path.exists(list_fqn):
        os.unlink(list_fqn)
    listobs(vis=fqn, listfile=list_fqn)


# cardinality == 1 obs per MOUS ID
# def test_gen():
#     # from astroquery.alma import Alma
#     # from astropy.table import Table
#     # db_table = Alma().query(payload={'project_code': '2016.1.00010.S'})
#     # # db_table.write('{}/alma_query.xml'.format(TEST_DATA_DIR),
#     # # format='votable')
#     # db_table.write('{}/alma_query2.xml'.format(TEST_DATA_DIR), format='html')
#
#     # FYI - this has values
#     # import logging
#     # logging.error(db_table['Release date'])
#     # assert False
#
#     import logging
#     from caom2.diff import get_differences
#     from caom2pipe import manage_composable as mc
#     from almaca2caom2 import main_app as ma
#
#     from astropy.table import Table
#     db_fqn = '{}/alma_query.html'.format(ALMA_QUERY_DIR)
#     db_content = Table.read(db_fqn, format='html')
#
#     obs = None
#     entries = glob.glob('/data/calibrated/*ms.split.cal')
#     for entry in entries:
#         fqn = f'{entry}/md.pk'
#         if not os.path.exists(fqn):
#             logging.error('missing md for {}'.format(fqn))
#             continue
#         logging.error('reading from {}'.format(fqn))
#         temp = mc.read_from_file(fqn)
#         assert temp is not None, 'expected result'
#         result = {'spectral_windows': []}
#         for line in temp:
#             temp2 = line.split(',', 1)
#             if temp2[0].strip() == 'spectral_windows':
#                 x = temp2[1].strip().split(',')
#                 count = 0
#                 while count < len(x):
#                     y = (mc.to_float(x[count].replace('(', '').replace('[', '').replace(']', '').replace(')', '')),
#                          mc.to_float(x[count+1].replace('(', '').replace('[', '').replace(']', '').replace(')', '')))
#                     count += 2
#                     result[temp2[0]].append(y)
#             elif type(temp2[1]) is str:
#                 result[temp2[0]] = temp2[1].strip()
#             else:
#                 result[temp2[0]] = temp2[1]
#         try:
#             obs = ma.build_observation(result, db_content, obs, fqn)
#         except Exception as e:
#             import traceback
#             logging.error(traceback.format_exc())
#             logging.error(fqn)
#             assert False
#
#     obs_fqn = f'{ALMA_QUERY_DIR}/{obs.observation_id}.actual.xml'
#     logging.error(f'obs type {type(obs)}')
#     mc.write_obs_to_file(obs, obs_fqn)
#     expected_fqn = '{}/expected_{}.xml'.format(ALMA_QUERY_DIR,
#                                                obs.observation_id)
#     if os.path.exists(expected_fqn):
#         expected = mc.read_obs_from_file(expected_fqn)
#         result = get_differences(expected, obs, 'Observation')
#         if result:
#             msg = 'Differences found {}\n{}'.format(
#                 obs.observation_id, '\n'.join([r for r in result]))
#             logging.error(msg)
#             errors_found = True
#             # assert False, obs_fqn
#     else:
#         raise AssertionError(
#             'Unexpected Observation {}'.format(obs.observation_id))
#
#     # assert False


# For a different cardinality handling - SGo - 02-03-20
def test_gen():
    # from astroquery.alma import Alma
    # from astropy.table import Table
    # db_table = Alma().query(payload={'project_code': '2016.1.00010.S'})
    # # db_table.write('{}/alma_query.xml'.format(TEST_DATA_DIR),
    # # format='votable')
    # db_table.write('{}/alma_query2.xml'.format(TEST_DATA_DIR), format='html')

    # FYI - this has values
    # import logging
    # logging.error(db_table['Release date'])
    # assert False

    import logging
    from caom2.diff import get_differences
    from caom2pipe import manage_composable as mc
    from almaca2caom2 import main_app as ma

    from astropy.table import Table
    db_fqn = '{}/alma_query.html'.format(ALMA_QUERY_DIR)
    db_content = Table.read(db_fqn, format='html')

    obs_lookup = {}
    obs = None
    entries = glob.glob('/data/calibrated/*ms.split.cal')
    for entry in entries:
        fqn = f'{entry}/md.pk'
        if not os.path.exists(fqn):
            logging.error('missing md for {}'.format(fqn))
            continue

        logging.error(entry)
        try:
            almaca_name = ma.AlmacaName(os.path.basename(entry))
        except mc.CadcException:
            # skip the input MS's
            continue

        if obs is not None and almaca_name.obs_id != obs.observation_id:
            # logging.error('storing {}'.format(obs.observation_id))
            obs_lookup[obs.observation_id] = obs
            obs = obs_lookup.get(almaca_name.obs_id)

        # temp = mc.read_from_file(fqn)
        # assert temp is not None, 'expected result'
        # result = {'spectral_windows': []}
        # for line in temp:
        #     temp2 = line.split(',', 1)
        #     if temp2[0].strip() == 'spectral_windows':
        #         x = temp2[1].strip().split(',')
        #         count = 0
        #         while count < len(x):
        #             y = (mc.to_float(x[count].replace('(', '').replace('[', '').replace(']', '').replace(')', '')),
        #                  mc.to_float(x[count+1].replace('(', '').replace('[', '').replace(']', '').replace(')', '')))
        #             count += 2
        #             result[temp2[0]].append(y)
        #     elif type(temp2[1]) is str:
        #         result[temp2[0]] = temp2[1].strip()
        #     else:
        #         result[temp2[0]] = temp2[1]
        try:
            obs = ma.build_observation(db_content, obs, fqn)
        except Exception as e:
            import traceback
            logging.error(traceback.format_exc())
            logging.error(fqn)
            assert False

    errors_found = False
    for obs in obs_lookup:
        obs_fqn = f'{ALMA_QUERY_DIR}/{obs}.actual.xml'
        logging.error(f'checking fqn {obs_fqn}')
        actual = obs_lookup.get(obs)
        mc.write_obs_to_file(actual, obs_fqn)
        expected_fqn = f'{ALMA_QUERY_DIR}/{actual.observation_id}.expected.xml'
        if os.path.exists(expected_fqn):
            expected = mc.read_obs_from_file(expected_fqn)
            result = get_differences(expected, actual, 'Observation')
            if result:
                msg_str = '\n'.join([r for r in result])
                msg = f'Differences found {actual.observation_id}\n{msg_str}'
                logging.error(msg)
                errors_found = True
                # assert False, obs_fqn
        else:
            msg = f'Could not find {expected_fqn}'
            logging.error(msg)
            raise AssertionError(msg)

    if errors_found:
        raise AssertionError('Final failure.')
    # assert False


def get_info(filename, provenance):
    import logging
    info = {}
    msmd.open(filename)
    info['band'] = int(msmd.namesforspws(0)[0].split("#")[1].split("_")[2])
    temp = msmd.fieldnames()
    if len(temp) > 1:
        # HK 13-08-2019
        # Given that, it may be most practical to simply check that all of
        # the 'fieldnames' entries are the same, and if they are, continue
        # as usual with the metadata harvesting.  If it were straightforward
        # to do, it might also be worth generating some kind of flag / note
        # / separate file about the fact that the field is given 3 times
        # rather than 1, in case this is something that we need to bring up
        # with the ALMA archive team at a later date.
        fields = list(set(temp))
        logging.warning('Repeated field names for MS {}'.format(filename))
        if len(fields) > 1:
            raise ValueError("Expected ms with just one target")
        field = fields[0]
    else:
        field = temp[0]
    if msmd.nobservations() > 1:
        raise ValueError("Exepected ms with just one observation")
    info['field'] = field
    logging.error('field is {}'.format(field))
    spws = msmd.spwsforfield(field)
    spectral_windows = []
    for idx in spws:
        spectral_windows.append((min(msmd.chanfreqs(idx, 'Hz')),
                                 max(msmd.chanfreqs(idx, 'Hz'))))
    info['spectral_windows'] = spectral_windows
    dates = msmd.timerangeforobs(0)
    info['start_date'] = dates['begin']['m0']['value']
    info['end_date'] = dates['end']['m0']['value']
    info['project'] = msmd.projects()[0]
    scans = msmd.scansforfield(field)
    itime = 0
    for scan in scans:
        itime += len(msmd.timesforscan(scan)) * msmd.exposuretime(scan)['value']
    info['itime'] = itime
    info['effexposuretime'] = msmd.effexposuretime()['value']

    # HK 19-08-19
    # energy: sample size - if this is the number of spectral resolution
    # elements, we could add this!  Adding up all of the elements in the
    # output for msmd.chanwidths() (run over each of the 4 spectral
    # windows) would give the answer here.

    # HK 19-08-19
    # energy: resolution - I don't see this listed in the caom2 model parameters
    # (http://www.opencadc.org/caom2/#Energy), but would this represent the
    # frequency resolution?  If so, it could be added in, with the
    # complication that this will be variable for the different spectral
    # windows, and hence one merged range of frequencies might have multiple
    # frequency resolutions.  The msmd.chanwidths function (run per spectral
    # window) will give this information in frequency space

    # PD - 05-09-19
    # resolution may be the same as sampleSize or it could be larger (photons
    # scattered into neighbouring elements); for most of our data time
    # resolution = sampleSize... might not be true for real time series
    # data... note that energy axis is the resolvingPower (unitless)
    #
    # sampleSize is the size of one element (aka pixel):
    #   - wavelength in m for energy
    #
    # HK 05-09-19
    # For the energy dimension, it sounds like we want the following mappings:
    #
    # - sampleSize --> msmd.chanwidths
    #
    # - resolution --> msmd.chanres
    #
    # (sometimes observations are set up with the expectation of
    # binning/smoothing, which would be captured in chanres but not
    # chanwidths.  For the observation I've been looking at, the two values
    # differ by a factor of 2, i.e., binning every 2 channels is expected)

    # HK 03-12-19
    # For the checks/discussion below, I am assuming the caom2 Energy listings
    # are in metres (except ResolvingPower which is unitless). As we discussed
    # earlier, SampleSize = msmd.chanwidth, Resolution=msmd.chanres (barring
    # unit conversions)
    #
    # I'm assuming that ResolvingPower = mean[frequency_Hz] / chanres_Hz
    #
    # Conversion of other parameters: since ResolvingPower should be the same
    # for either frequency or wavelength units, we get:
    #
    # chanwidth_metres = (chanwidth_Hz/ mean[frequency_Hz]) * mean[wavelength_m]
    #
    # and similarly for chanres_metres

    # assuming 0th index, assuming the widths are the same for all
    chan_widths = msmd.chanwidths(0, 'Hz')
    temp = list(set(chan_widths))
    if len(temp) > 1:
        logging.warning('Found multiple values for chan_widths')
    else:
        info['energy_sample_size'] = chan_widths[0]

    # assuming 0th index, assuming the res is the same for all
    chan_res = msmd.chanres(0, 'Hz', asvel=False)
    temp = list(set(chan_res))
    if len(temp) > 1:
        logging.warning('Found multiple values for chan_res')
    else:
        info['energy_resolution'] = temp[0]

    # for idx in spws:
    #     logging.error('idx is {} len is {}'.format(idx, len(msmd.chanwidths(idx))))
    #     logging.error(msmd.chanwidths(idx))
    #
    #     # info['energy_resolution'] = msmd.chanres(idx)
    #     logging.error('chanres len is {}'.format(len(msmd.chanres(idx))))
    #     logging.error(msmd.chanres(idx))


    # PD - 05-09-19
    # sampleSize is the size of one element (aka pixel):
    #   - time in days for time axis
    #
    # HK 05-09-19
    # For the time dimension:
    #
    # - sampleSize --> # of elements in msmd.timesforfield
    # got through all the elements identified by msmd.fieldsfortimes, not just
    # the 0th one
    fields_for_times = msmd.fieldsfortimes()
    sample_size = 0
    for field in fields_for_times:
        sample_size += len(msmd.timesforfield(field))
    info['time_sample_size'] = sample_size
    info['provenance'] = provenance

    msmd.done()
    # testing below
    # print('bandwidths')
    # print(msmd.bandwidths())
    # print('chanfreqs')
    # for ii in spws:
    #     print(msmd.chanfreqs(ii))
    # print(dir(msmd))
    # try:
    #     for ii in dir(msmd):
    #         # the named entries just need parameters
    #         if not ii.startswith('_') and ii not in [  # 'antennadiameter',
    #                                                  'antennaoffset',
    #                                                  'antennasforscan',
    #                                                  'baseband',
    #                                                  'chaneffbws',
    #                                                  'chanfreqs',
    #                                                  'chanavgspws',
    #                                                  'chanres',
    #                                                  'chanwidths',
    #                                                  'corrprodsforpol',
    #                                                  'corrtypesforpol',
    #                                                  'close',
    #                                                  'datadescids',
    #                                                  'done',
    #                                                  # 'effexposuretime',
    #                                                  'exposuretime',
    #                                                  # 'fdmspws',
    #                                                  # 'fieldnames',
    #                                                  # 'fieldsforintent',
    #                                                  # 'fieldsforname',
    #                                                  'fieldsforscan',
    #                                                  # 'fieldsforscans',
    #                                                  # 'fieldsforsource',
    #                                                  # 'fieldsforsources',
    #                                                  'fieldsforspw',
    #                                                  # 'fieldsfortimes',
    #                                                  # 'intents',
    #                                                  'intentsforfield',
    #                                                  'intentsforscan',
    #                                                  'intentsforspw',
    #                                                  'meanfreq',
    #                                                  # 'name',
    #                                                  # 'namesforfields',
    #                                                  # 'namesforspws',
    #                                                  # 'nantennas',
    #                                                  # 'narrays',
    #                                                  # 'nbaselines',
    #                                                  'nchan',
    #                                                  # 'ncorrforpol',
    #                                                  # 'nfields',
    #                                                  # 'nobservations',
    #                                                  # 'nrows',
    #                                                  # 'nscans',
    #                                                  # 'nsources',
    #                                                  # 'nspw',
    #                                                  # 'nstates',
    #                                                  # 'observatorynames',
    #                                                  # 'observatoryposition',
    #                                                  # 'observers',
    #                                                  'open',
    #                                                  # 'phasecenter',
    #                                                  # 'pointingdirection',
    #                                                  # 'polidfordatadesc',
    #                                                  # 'projects',
    #                                                  # 'propermotions',
    #                                                  'refdir',
    #                                                  'reffreq',
    #                                                  # 'restfreqs',
    #                                                  # 'scannumbers',
    #                                                  # 'scansforfields',
    #                                                  'scansforfield',
    #                                                  'scansforintent',
    #                                                  'scansforspw',
    #                                                  'scansforstate',
    #                                                  'schedule',
    #                                                  'sideband',
    #                                                  'sourceidforfield',
    #                                                  'spwsforfield',
    #                                                  'spwsforscan',
    #                                                  'statesforscan',
    #                                                  'this',
    #                                                  'timerangeforobs',
    #                                                  'timesforscan',
    #                                                  'timesforscans'
    #                                                  ]:
    #             print(ii)
    #             print(getattr(msmd, ii)())
    # except Exception as e:
    #     import traceback
    #     print(traceback.format_exc())
    return info


if __name__ == '__main__':
    run_gen()
    # run_gen_specific()
