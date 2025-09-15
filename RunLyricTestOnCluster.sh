## To run lyric test on cluster
git clone https://github.com/guigolab/LyRic lyric_test
cd lyric_test

module load Python/3.13.5-GCCcore-14.3.0
python3 -m venv snakemake8_lyrictest
source snakemake8_lyrictest/bin/activate
cd snakemake8_lyrictest/bin/
./pip install snakemake==8 pandas pkg_resources

cd ../../
snakemake --cores 1