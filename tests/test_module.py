# Test the function
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

# Importing the functions from molinterface
from molinterface import (
    find_smiles_patterns,
    find_rings,
    format_rings,
    get_compound_info,
    highlight_chemical_groups,
    find_chemical_groups,
    open_url
)

def test_find_smiles_patterns():
    smiles = 'CCO'
    smarts_patterns = {'[CX3H2][OX2H]': 'Primary Alcohol'}
    result = find_smiles_patterns(smiles, smarts_patterns)
    assert 'Primary Alcohol' in result

def test_find_rings():
    smiles = 'C1=CC=CC=C1'  # Benzene
    result = find_rings(smiles)
    assert result['non fused aromatic rings']['6 membered ring'] == 1

def test_format_rings():
    rings_dict = {
        'non fused aromatic rings': {'6 membered ring': 1},
        'non fused non aromatic rings': {},
        'fused aromatic rings': {},
        'fused non aromatic rings': {},
        'fused mixed rings': {}
    }
    result = format_rings(rings_dict)
    expected = 'Non fused aromatic rings:  6 membered ring: 1 | '
    assert result == expected

def test_get_compound_info():
    molecule_name = 'ethanol'
    smarts_patterns = {'[CX3H2][OX2H]': 'Primary Alcohol'}
    result = get_compound_info(molecule_name, smarts_patterns)
    assert 'iupac_name' in result
    assert result['iupac_name'] == 'ethanol'

def test_highlight_chemical_groups():
    smiles = 'CCO'
    smarts_patterns = {'[CX3H2][OX2H]': 'Primary Alcohol'}
    img_data = highlight_chemical_groups(smiles, smarts_patterns)
    assert isinstance(img_data, bytes)  # Check if the result is in bytes (Cairo image data)

def test_find_chemical_groups():
    smiles = 'CCO'
    smarts_patterns = {'[CX3H2][OX2H]': 'Primary Alcohol'}
    result = find_chemical_groups(smiles, smarts_patterns)
    assert 'Primary Alcohol' in result

def test_open_url(monkeypatch):
    import webbrowser
    url = 'http://example.com'

    def mock_open_new(url):
        return True

    monkeypatch.setattr(webbrowser, 'open_new', mock_open_new)
    open_url(url)
    assert webbrowser.open_new.called_once_with(url)

if __name__ == '__main__':
    pytest.main()

