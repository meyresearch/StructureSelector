
import StructureSelector.prepper as prepper
import os


def test_download_pdb():
    p = prepper.Prepper()
    p.download_pdbs(['1FKN'], '/tmp/')
    assert (os.path.isfile('/tmp/1FKN.pdb'))

def test_download_pdb_httpError():
    print('Not implemented yet')