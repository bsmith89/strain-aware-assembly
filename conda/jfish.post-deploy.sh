# NOTE: While this script is correct in principle, it seems to fail
# on the linking because $CONDA_PREFIX is not set correctly by Snakemake.
# history > jfish.post-deploy.log
# env >> jfish.post-deploy.log
# echo $CONDA_PREFIX >> jfish.post-deploy.log
# for basename in dna_jellyfish-0.0.1.dist-info dna_jellyfish.py _dna_jellyfish.cpython-39-x86_64-linux-gnu.so
# do
#     ln -s /opt/conda/lib/python3.9/site-packages/$basename $CONDA_PREFIX/lib/python3.9/site-packages/$basename
# done
# python -c "import dna_jellyfish" >> jfish.post-deploy.log 2>&1
