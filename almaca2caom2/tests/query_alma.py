from astroquery.alma import Alma
from astropy.table import Table

# db_table = Alma().query(payload={'project_code': '2016.1.00010.S'})
db_table = Alma().query(payload={'project_code': '2016.1.00010.S'},
                        result_view='project')
# db_table.write('{}/alma_query.xml'.format(TEST_DATA_DIR),
# format='votable')
db_table.write('{}/alma_query3.html'.format('./'), format='html', overwrite=True)
