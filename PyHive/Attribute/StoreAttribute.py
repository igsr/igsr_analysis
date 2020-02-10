import eHive

from ReseqTrackDB import Attribute
from ReseqTrackDB import ReseqTrackDB


class StoreAttribute(eHive.BaseRunnable):
    """Store an Attribute or a list of Attributes in the DB"""

    def fetch_input(self):
        hostname = self.param_required('hostname')
        username = self.param_required('username')
        db = self.param_required('db')
        port = self.param_required('port')
        pwd = self.param_required('pwd')

        reseqdb = ReseqTrackDB(host=hostname, user=username, port=port, pwd=pwd, db=db)

        self.param('reseqdb', reseqdb)

    def run(self):

        attrb = self.param_required('attrb')
        reseqdb = self.param('reseqdb')

        attO = Attribute(attrb['table_name'], attrb['other_id'], attrb['name'], attrb['value'])

        self.warning('Attribute with name: {0}'.format(attO.name))
        if self.param_required('store_attributes') == 'True':
            attO.store(reseqdb)
            self.warning('Attribute with name: {0} was stored in DB'.format(attO.name))

    def write_output(self):
        self.warning('Work is done!')
