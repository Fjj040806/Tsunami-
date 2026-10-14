"""Reproduce the full study; synthetic-only mode requires no external data."""
from pathlib import Path
import argparse,subprocess,sys
ROOT=Path(__file__).resolve().parents[1]
def run(script,*args):
    print('\nRunning',script,*args,flush=True)
    subprocess.run([sys.executable,str(ROOT/'scripts'/script),*args],cwd=ROOT,check=True)
def main():
    p=argparse.ArgumentParser();p.add_argument('--synthetic-only',action='store_true');a=p.parse_args()
    run('run_verification.py');run('run_experiments.py');run('run_v02.py');run('run_v03.py')
    for part in ['step','shelf']:run('run_final_studies.py','--section',part)
    for ratio in ['0.5','0.25']:run('run_weak_refinement.py','--ratio',ratio)
    if not a.synthetic_only:
        run('run_final_studies.py','--section','lab');run('build_final_figures.py')
    else:
        print('Laboratory comparison and final report-table builder skipped; bundled final report figures remain available.')
    run('build_visualizations.py')
    print('Reproduction completed. To rebuild the PDF, run pdflatex twice from paper/.')
if __name__=='__main__':main()
