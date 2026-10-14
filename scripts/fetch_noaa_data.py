"""Download the NOAA-hosted USACE archive; fail if its pinned bytes change."""
from pathlib import Path
import urllib.request,hashlib,json
ROOT=Path(__file__).resolve().parents[1]
URL='https://nctr.pmel.noaa.gov/benchmark/Solitary_wave/tsunami3_runup.zip'
def main():
    p=ROOT/'data/raw/noaa_composite.zip';p.parent.mkdir(parents=True,exist_ok=True)
    expected=json.loads((ROOT/'data/noaa_provenance.json').read_text())['sha256']
    data=urllib.request.urlopen(URL,timeout=60).read()
    if hashlib.sha256(data).hexdigest()!=expected:raise RuntimeError('Dataset hash changed; audit before use')
    p.write_bytes(data);print('NOAA archive downloaded and checksum verified.')
if __name__=='__main__':main()
