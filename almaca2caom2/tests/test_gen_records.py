import os
import math

TEST_DATA_DIR = '/usr/src/app/almaca2caom2/almaca2caom2/tests/data'


def run_gen():
    import logging
    for root, dirs, files in os.walk(TEST_DATA_DIR):
        for dir_name in dirs:
            if 'ms.split.cal' in dir_name:
                fqn = '{}/{}'.format(root, dir_name)
                try:
                    logging.error('Processing {}'.format(fqn))
                    result = get_info(fqn)
                    # result = try_vishead(fqn)
                    logging.error(result)
                    assert result is not None, 'expected result'
                    pk_file = '{}/md.pk'.format(TEST_DATA_DIR)
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
    # import logging
    # taskname='vishead'
    # default()
    # vis=fqn
    # mode='get'
    # hdkey='field'
    # hdindex='2'
    # hdvalue=vishead()
    # logging.error(hdvalue[0])
    listobs(vis=fqn, listfile='./abc.txt')


def test_gen():
    # from astroquery.alma import Alma
    # from astropy.table import Table
    # db_table = Alma().query(payload={'project_code': '2016.1.00010.S'})
    # db_table.write('{}/alma_query.xml'.format(TEST_DATA_DIR), format='votable')
    # assert False

    import logging
    from caom2.diff import get_differences
    from caom2pipe import manage_composable as mc
    from almaca2caom2 import main_app as ma

    from astropy.table import Table
    db_fqn = '{}/alma_query.xml'.format(TEST_DATA_DIR)
    db_content = Table.read(db_fqn, format='votable')

    sci_obs = None
    cal_obs = None
    for ii in os.listdir(TEST_DATA_DIR):
        if ii.endswith('.pk'):
            if 'cal' in ii:
                obs = cal_obs
            else:
                obs = sci_obs
            fqn = '{}/{}'.format(TEST_DATA_DIR, ii)
            logging.error(fqn)
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
            try:
                obs = ma.build_observation(result, db_content, obs, fqn)
            except Exception as e:
                import traceback
                logging.error(traceback.format_exc())
                assert False
            if 'cal' in ii:
                cal_obs = obs
            else:
                sci_obs = obs

    for obs in [cal_obs, sci_obs]:
        obs_fqn = '{}/actual_{}.xml'.format(TEST_DATA_DIR, obs.observation_id)
        mc.write_obs_to_file(obs, obs_fqn)
        expected_fqn = '{}/expected_{}.xml'.format(TEST_DATA_DIR, obs.observation_id)
        expected = mc.read_obs_from_file(expected_fqn)
        result = get_differences(expected, obs, 'Observation')
        if result:
            logging.error('From {}'.format(expected_fqn))
            msg = 'Differences found\n{}'.format('\n'.join([r for r in result]))
            raise AssertionError(msg)
    # assert False


def get_info(filename):
    import logging
    info = {}
    msmd.open(filename)
    logging.error(msmd.namesforspws())
    info['band'] = int(msmd.namesforspws(0)[0].split("#")[1].split("_")[2])
    position = msmd.phasecenter()
    info['ra'] = math.degrees(position['m0']['value'])
    info['dec'] = math.degrees(position['m1']['value'])
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
    else:
        fields = temp
    if msmd.nobservations() > 1:
        raise ValueError("Exepected ms with just one observation")
    info['field'] = fields[0]
    spws = msmd.spwsforfield(fields[0])
    spectral_windows = []
    for idx in spws:
        spectral_windows.append((min(msmd.chanfreqs(idx)), max(msmd.chanfreqs(idx))))
    info['spectral_windows'] = spectral_windows
    dates = msmd.timerangeforobs(0)
    info['start_date'] = dates['begin']['m0']['value']
    info['end_date'] = dates['end']['m0']['value']
    info['project'] = msmd.projects()[0]
    scans = msmd.scansforfield(info['field'])
    itime = 0
    itime_len = 0
    for scan in scans:
        itime += len(msmd.timesforscan(scan)) * msmd.exposuretime(scan)['value']
        itime_len += len(msmd.timesforscan(scan))
    info['itime'] = itime
    info['resolution'] = itime_len
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
    # window) will give this information in frequency space - this is still
    # TODO until there's a resolution on the whether or not to split on
    # Energy Band, because that will provide an answer as to 'which value to
    # use'
    sample_size = 0
    energy_resolution = 0
    for ii in spws:
        logging.error('idx is {} value is {}'.format(ii, temp[0]))
        kwargs = {'spw': ii}
        temp = msmd.chanwidths(**kwargs)
        sample_size += len(temp)
    info['sample_size'] = sample_size
    info['energy_resolution'] = energy_resolution

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
    #                                                  'timesforfield',
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
