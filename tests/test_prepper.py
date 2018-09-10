
import StructureSelector.prepper as prepper
import os
import requests
import pytest


def test_download_pdb():
    p = prepper.Prepper()
    p.download_pdbs(['1FKN'], '/tmp/')
    assert (os.path.isfile('/tmp/1FKN.pdb'))

def test_download_pdb_httpError():
    p = prepper.Prepper()
    p.download_pdbs(['XXXX'], '/tmp/')
    assert(pytest.raises(requests.exceptions.HTTPError))

def test_clean_pdb():
    pass

def test_split_chains():
    pass

def test_pdb_to_fasta():
    pass