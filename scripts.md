## RNA-Seq (k=23)

```bash
for N in {1..10}; do
    bsub -J "subset_${N}" \
         -oo ~/metagenome/data/mantis/subsets/mtg/canonical/graph_canonical_subset_${N}k.lsf \
         -W 8:00 \
         -n 15 -R "rusage[mem=22000] span[hosts=1]" \
        "cat ~/metagenome/data/mantis/subsets/${N}k.txt \
            | /usr/bin/time -v ~/projects/projects2014-metagenome/metagraph/build_test/metagraph build -v \
                -k 23 \
                --canonical \
                --parallel 30 \
                --mem-cap-gb 300 \
                --disk-cap-gb 10000 \
                --disk-swap ~/metagenome/data/mantis/subsets/mtg/canonical/ \
                -o ~/metagenome/data/mantis/subsets/mtg/canonical/graph_canonical_subset_${N}k \
                2>&1"; \
done


for N in {1..10}; do
    bsub -J "subset_transform_${N}" \
         -oo ~/metagenome/data/mantis/subsets/mtg/canonical/graph_canonical_subset_${N}k_transform.lsf \
         -w "subset_${N}" \
         -W 8:00 \
         -n 15 -R "rusage[mem=22000] span[hosts=1]" \
        "/usr/bin/time -v ~/projects/projects2014-metagenome/metagraph/build_test/metagraph transform -v \
                --to-fasta \
                --parallel 30 \
                --primary-kmers \
                -o ~/metagenome/data/mantis/subsets/mtg/canonical/graph_primary_subset_${N}k \
                ~/metagenome/data/mantis/subsets/mtg/canonical/graph_canonical_subset_${N}k.dbg \
                2>&1"; \
done


for N in {1..10}; do
    bsub -J "primary_${N}" \
         -oo ~/metagenome/data/mantis/subsets/mtg/logs/graph_primary_subset_${N}k_build.lsf \
         -w "subset_transform_${N}" \
         -W 8:00 \
         -n 15 -R "rusage[mem=22000] span[hosts=1]" \
        "/usr/bin/time -v ~/projects/projects2014-metagenome/metagraph/build_test/metagraph build -v \
                -k 23 \
                --parallel 30 \
                --mem-cap-gb 300 \
                --disk-cap-gb 10000 \
                --disk-swap ~/metagenome/data/mantis/subsets/mtg/ \
                -o ~/metagenome/data/mantis/subsets/mtg/graph_primary_subset_${N}k \
                ~/metagenome/data/mantis/subsets/mtg/canonical/graph_primary_subset_${N}k.fasta.gz \
                2>&1"; \
done


for N in {1..10}; do
    mkdir ~/metagenome/data/mantis/subsets/mtg/cols_${N}k;
    bsub -J "annotate_${N}" \
         -oo ~/metagenome/data/mantis/subsets/mtg/logs/graph_primary_subset_${N}k_annotate.lsf \
         -w "primary_${N}" \
         -W 8:00 \
         -n 15 -R "rusage[mem=22000] span[hosts=1]" \
        "cat ~/metagenome/data/mantis/subsets/${N}k.txt \
            | /usr/bin/time -v ~/projects/projects2014-metagenome/metagraph/build_test/metagraph annotate -v \
                --parallel 30 \
                --canonical \
                --anno-filename \
                --separately \
                -o ~/metagenome/data/mantis/subsets/mtg/cols_${N}k \
                -i ~/metagenome/data/mantis/subsets/mtg/graph_primary_subset_${N}k.dbg \
                2>&1"; \
done


DIR=~/metagenome/data/row_diff/subsets/rnaseq/nobackup/mtg
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_master/metagraph
for DEPTH in {100,75,50,25,10}; do
    THREADS=72
    for N in {10,}; do
        NAME=rd_${DEPTH}d_${N}k_${THREADS}t;
        RD_DIR=${DIR}/${NAME};
        rm -r ${RD_DIR};
        mkdir ${RD_DIR};
        ln -s ${DIR}/graph_primary_subset_${N}k.dbg \
              ${RD_DIR}/graph.dbg;
        bsub -J "${NAME}" \
             -oo ${DIR}/logs/${NAME}.lsf \
             -W 10:00 \
             -n 36 -R "rusage[mem=9500] span[hosts=1]" \
            "find ${DIR}/cols_${N}k -name \"*.column.annodbg\" \
                | /usr/bin/time -v $METAGRAPH transform_anno -v \
                    --anno-type row_diff \
                    --max-path-length $DEPTH \
                    --mem-cap-gb 300 \
                    --parallel $THREADS \
                    --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
                    -o ${RD_DIR}/out \
                    -i ${RD_DIR}/graph.dbg \
                    2>&1 | tee ${DIR}/logs/${NAME}.log";
        bsub -J "${NAME}_opt" \
             -w "${NAME}" \
             -oo ${DIR}/logs/${NAME}_opt.lsf \
             -W 10:00 \
             -n 36 -R "rusage[mem=9500] span[hosts=1]" \
            "find ${DIR}/cols_${N}k -name \"*.column.annodbg\" \
                | /usr/bin/time -v $METAGRAPH transform_anno -v \
                    --anno-type row_diff \
                    --max-path-length $DEPTH \
                    --mem-cap-gb 600 \
                    --parallel $THREADS \
                    --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
                    --optimize \
                    -o ${RD_DIR}/out \
                    -i ${RD_DIR}/graph.dbg \
                    2>&1 | tee ${DIR}/logs/${NAME}_opt.log"; \
        bsub -J "${NAME}_multi_brwt" \
             -w "${NAME}_opt" \
             -oo ${DIR}/logs/${NAME}_multi_brwt.lsf \
             -W 10:00 \
             -n 36 -R "rusage[mem=9500] span[hosts=1]" \
            "find ${RD_DIR} -name \"*.row_diff.annodbg\" \
                | /usr/bin/time -v $METAGRAPH transform_anno -v \
                    --anno-type row_diff_brwt \
                    --greedy \
                    --parallel $THREADS \
                    --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
                    -o ${DIR}/multi_brwt/${NAME} \
                    -i ${RD_DIR}/graph.dbg \
                    2>&1 | tee ${DIR}/logs/${NAME}_multi_brwt.log"; \
        bsub -J "${NAME}_multi_brwt_relax" \
             -w "${NAME}_multi_brwt" \
             -oo ${DIR}/logs/${NAME}_multi_brwt_relax.lsf \
             -W 10:00 \
             -n 36 -R "rusage[mem=9500] span[hosts=1]" \
            "/usr/bin/time -v $METAGRAPH relax_brwt -v \
                    --relax-arity 20 \
                    --parallel $THREADS \
                    --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
                    -o ${DIR}/multi_brwt/${NAME}_relaxed \
                    ${DIR}/multi_brwt/${NAME}.row_diff_brwt.annodbg \
                    2>&1 | tee ${DIR}/logs/${NAME}_multi_brwt_relax.log"; \
        bsub -J "${NAME}_row_sparse" \
             -w "${NAME}_opt" \
             -oo ${DIR}/logs/${NAME}_row_sparse.lsf \
             -W 10:00 \
             -n 36 -R "rusage[mem=9500] span[hosts=1]" \
            "find ${RD_DIR} -name \"*.row_diff.annodbg\" \
                | /usr/bin/time -v $METAGRAPH transform_anno -v \
                    --anno-type row_diff_sparse \
                    --parallel $THREADS \
                    --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
                    -o ${DIR}/row_sparse/${NAME} \
                    -i ${RD_DIR}/graph.dbg \
                    2>&1 | tee ${DIR}/logs/${NAME}_row_sparse.log"; \
    done
done


for N in {1..10}; do
    bsub -J "rb_brwt_${N}" \
         -oo ~/metagenome/data/mantis/subsets/mtg/logs/graph_primary_subset_${N}k_rb_brwt.lsf \
         -W 48:00 \
         -n 15 -R "rusage[mem=41000] span[hosts=1]" \
        "find ~/metagenome/data/mantis/subsets/mtg/cols_${N}k -name \"*.column.annodbg\" \
            | /usr/bin/time -v ~/projects/projects2014-metagenome/metagraph/build_test/metagraph transform_anno -v \
                --anno-type rb_brwt \
                --greedy --relax-arity 20 \
                --parallel 30 \
                -o ~/metagenome/data/mantis/subsets/mtg/rb_brwt_${N}k \
                2>&1 | tee ~/metagenome/data/mantis/subsets/mtg/logs/graph_primary_subset_${N}k_rb_brwt.log"; \
done


for N in {1..10}; do
    bsub -J "brwt_${N}" \
         -w "annotate_${N}" \
         -oo ~/metagenome/data/mantis/subsets/mtg/logs/graph_primary_subset_${N}k_brwt.lsf \
         -W 48:00 \
         -n 15 -R "rusage[mem=21000] span[hosts=1]" \
        "find ~/metagenome/data/mantis/subsets/mtg/cols_${N}k -name \"*.column.annodbg\" \
            | /usr/bin/time -v ~/projects/projects2014-metagenome/metagraph/build_test/metagraph transform_anno -v \
                --anno-type brwt \
                --greedy \
                --parallel 30 \
                -o ~/metagenome/data/mantis/subsets/mtg/brwt_${N}k \
                2>&1 | tee ~/metagenome/data/mantis/subsets/mtg/logs/graph_primary_subset_${N}k_brwt.log"; \
    bsub -J "brwt_${N}_relax" \
         -w "brwt_${N}" \
         -oo ~/metagenome/data/mantis/subsets/mtg/logs/graph_primary_subset_${N}k_brwt_relax.lsf \
         -W 48:00 \
         -n 15 -R "rusage[mem=21000] span[hosts=1]" \
        "/usr/bin/time -v ~/projects/projects2014-metagenome/metagraph/build_test/metagraph relax_brwt -v \
                --relax-arity 20 \
                --parallel 30 \
                -o ~/metagenome/data/mantis/subsets/mtg/brwt_relaxed_${N}k \
                ~/metagenome/data/mantis/subsets/mtg/brwt_${N}k.brwt.annodbg \
                2>&1 | tee ~/metagenome/data/mantis/subsets/mtg/logs/graph_primary_subset_${N}k_brwt_relax.log"; \
done

```


## RNA-Seq (k=31)

```bash

for N in {1..10}; do
    bsub -J "clean_subset_${N}" \
         -oo ~/metagenome/data/mantis/subsets/clean_rnaseq/mtg/canonical/clean_canonical_subset_${N}k.lsf \
         -W 8:00 \
         -n 15 -R "rusage[mem=22000] span[hosts=1]" \
        "cat ~/metagenome/data/mantis/subsets/clean_rnaseq/clean_${N}k.txt \
            | /usr/bin/time -v ~/projects/projects2014-metagenome/metagraph/build_test/metagraph build -v \
                -k 31 \
                --canonical \
                --parallel 30 \
                --mem-cap-gb 300 \
                --disk-cap-gb 10000 \
                --disk-swap ~/metagenome/data/mantis/subsets/clean_rnaseq/mtg/canonical/ \
                -o ~/metagenome/data/mantis/subsets/clean_rnaseq/mtg/canonical/clean_canonical_subset_${N}k \
                2>&1"; \
done


for N in {1..10}; do
    bsub -J "clean_transform_${N}" \
         -oo ~/metagenome/data/mantis/subsets/clean_rnaseq/mtg/canonical/clean_canonical_subset_${N}k_transform.lsf \
         -w "clean_subset_${N}" \
         -W 8:00 \
         -n 15 -R "rusage[mem=22000] span[hosts=1]" \
        "/usr/bin/time -v ~/projects/projects2014-metagenome/metagraph/build_test/metagraph transform -v \
                --to-fasta \
                --parallel 30 \
                --primary-kmers \
                -o ~/metagenome/data/mantis/subsets/clean_rnaseq/mtg/canonical/clean_primary_subset_${N}k \
                ~/metagenome/data/mantis/subsets/clean_rnaseq/mtg/canonical/clean_canonical_subset_${N}k.dbg \
                2>&1"; \
done


for N in {1..10}; do
    bsub -J "clean_primary_${N}" \
         -oo ~/metagenome/data/mantis/subsets/clean_rnaseq/mtg/logs/clean_primary_subset_${N}k_build.lsf \
         -w "clean_transform_${N}" \
         -W 8:00 \
         -n 15 -R "rusage[mem=22000] span[hosts=1]" \
        "/usr/bin/time -v ~/projects/projects2014-metagenome/metagraph/build_test/metagraph build -v \
                -k 31 \
                --parallel 30 \
                --mem-cap-gb 300 \
                --disk-cap-gb 10000 \
                --disk-swap ~/metagenome/data/mantis/subsets/clean_rnaseq/mtg/ \
                -o ~/metagenome/data/mantis/subsets/clean_rnaseq/mtg/clean_primary_subset_${N}k \
                ~/metagenome/data/mantis/subsets/clean_rnaseq/mtg/canonical/clean_primary_subset_${N}k.fasta.gz \
                2>&1"; \
done


for N in {1..10}; do
    mkdir ~/metagenome/data/mantis/subsets/clean_rnaseq/mtg/cols_${N}k;
    bsub -J "clean_annotate_${N}" \
         -oo ~/metagenome/data/mantis/subsets/clean_rnaseq/mtg/logs/clean_primary_subset_${N}k_annotate.lsf \
         -w "clean_primary_${N}" \
         -W 15:00 \
         -n 30 -R "rusage[mem=15000] span[hosts=1]" \
        "cat ~/metagenome/data/mantis/subsets/clean_rnaseq/clean_${N}k.txt \
            | /usr/bin/time -v ~/projects/projects2014-metagenome/metagraph/build_test/metagraph annotate -v \
                --parallel 60 \
                --canonical \
                --anno-filename \
                --separately \
                -o ~/metagenome/data/mantis/subsets/clean_rnaseq/mtg/cols_${N}k \
                -i ~/metagenome/data/mantis/subsets/clean_rnaseq/mtg/clean_primary_subset_${N}k.dbg \
                2>&1 | tee ~/metagenome/data/mantis/subsets/clean_rnaseq/mtg/logs/clean_primary_subset_${N}k_annotate.log"; \
done


DIR=~/metagenome/data/row_diff/subsets/clean_rnaseq/nobackup/mtg
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_master/metagraph
for DEPTH in {100,75,50,25,10}; do
    THREADS=72
    for N in {10,}; do
        NAME=clean_rd_${DEPTH}d_${N}k_${THREADS}t;
        RD_DIR=${DIR}/${NAME};
        rm -r ${RD_DIR};
        mkdir ${RD_DIR};
        ln -s ${DIR}/clean_primary_subset_${N}k.dbg \
              ${RD_DIR}/graph.dbg;
        bsub -J "${NAME}" \
             -oo ${DIR}/logs/${NAME}.lsf \
             -W 10:00 \
             -n 36 -R "rusage[mem=9500] span[hosts=1]" \
            "find ${DIR}/cols_${N}k -name \"*.column.annodbg\" \
                | /usr/bin/time -v $METAGRAPH transform_anno -v \
                    --anno-type row_diff \
                    --max-path-length $DEPTH \
                    --mem-cap-gb 300 \
                    --parallel $THREADS \
                    --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
                    -o ${RD_DIR}/out \
                    -i ${RD_DIR}/graph.dbg \
                    2>&1 | tee ${DIR}/logs/${NAME}.log";
        bsub -J "${NAME}_opt" \
             -w "${NAME}" \
             -oo ${DIR}/logs/${NAME}_opt.lsf \
             -W 10:00 \
             -n 36 -R "rusage[mem=9500] span[hosts=1]" \
            "find ${DIR}/cols_${N}k -name \"*.column.annodbg\" \
                | /usr/bin/time -v $METAGRAPH transform_anno -v \
                    --anno-type row_diff \
                    --max-path-length $DEPTH \
                    --mem-cap-gb 600 \
                    --parallel $THREADS \
                    --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
                    --optimize \
                    -o ${RD_DIR}/out \
                    -i ${RD_DIR}/graph.dbg \
                    2>&1 | tee ${DIR}/logs/${NAME}_opt.log"; \
        bsub -J "${NAME}_multi_brwt" \
             -w "${NAME}_opt" \
             -oo ${DIR}/logs/${NAME}_multi_brwt.lsf \
             -W 10:00 \
             -n 36 -R "rusage[mem=9500] span[hosts=1]" \
            "find ${RD_DIR} -name \"*.row_diff.annodbg\" \
                | /usr/bin/time -v $METAGRAPH transform_anno -v \
                    --anno-type row_diff_brwt \
                    --greedy \
                    --parallel $THREADS \
                    --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
                    -o ${DIR}/multi_brwt/${NAME} \
                    -i ${RD_DIR}/graph.dbg \
                    2>&1 | tee ${DIR}/logs/${NAME}_multi_brwt.log"; \
        bsub -J "${NAME}_multi_brwt_relax" \
             -w "${NAME}_multi_brwt" \
             -oo ${DIR}/logs/${NAME}_multi_brwt_relax.lsf \
             -W 10:00 \
             -n 36 -R "rusage[mem=9500] span[hosts=1]" \
            "/usr/bin/time -v $METAGRAPH relax_brwt -v \
                    --relax-arity 20 \
                    --parallel $THREADS \
                    --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
                    -o ${DIR}/multi_brwt/${NAME}_relaxed \
                    ${DIR}/multi_brwt/${NAME}.row_diff_brwt.annodbg \
                    2>&1 | tee ${DIR}/logs/${NAME}_multi_brwt_relax.log"; \
        bsub -J "${NAME}_row_sparse" \
             -w "${NAME}_opt" \
             -oo ${DIR}/logs/${NAME}_row_sparse.lsf \
             -W 10:00 \
             -n 36 -R "rusage[mem=9500] span[hosts=1]" \
            "find ${RD_DIR} -name \"*.row_diff.annodbg\" \
                | /usr/bin/time -v $METAGRAPH transform_anno -v \
                    --anno-type row_diff_sparse \
                    --parallel $THREADS \
                    --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
                    -o ${DIR}/row_sparse/${NAME} \
                    -i ${RD_DIR}/graph.dbg \
                    2>&1 | tee ${DIR}/logs/${NAME}_row_sparse.log"; \
    done
done


for N in {1..10}; do
    bsub -J "clean_rb_brwt_${N}" \
         -w "clean_annotate_${N}" \
         -oo ~/metagenome/data/mantis/subsets/clean_rnaseq/mtg/logs/clean_primary_subset_${N}k_rb_brwt.lsf \
         -W 48:00 \
         -n 15 -R "rusage[mem=41000] span[hosts=1]" \
        "find ~/metagenome/data/mantis/subsets/clean_rnaseq/mtg/cols_${N}k -name \"*.column.annodbg\" \
            | /usr/bin/time -v ~/projects/projects2014-metagenome/metagraph/build_test/metagraph transform_anno -v \
                --anno-type rb_brwt \
                --greedy --relax-arity 20 \
                --parallel 30 \
                -o ~/metagenome/data/mantis/subsets/clean_rnaseq/mtg/rb_brwt_${N}k \
                2>&1 | tee ~/metagenome/data/mantis/subsets/clean_rnaseq/mtg/logs/clean_primary_subset_${N}k_rb_brwt.log"; \
done


for N in {1..10}; do
    bsub -J "clean_brwt_${N}" \
         -w "clean_annotate_${N}" \
         -oo ~/metagenome/data/mantis/subsets/clean_rnaseq/mtg/logs/clean_primary_subset_${N}k_brwt.lsf \
         -W 48:00 \
         -n 15 -R "rusage[mem=41000] span[hosts=1]" \
        "find ~/metagenome/data/mantis/subsets/clean_rnaseq/mtg/cols_${N}k -name \"*.column.annodbg\" \
            | /usr/bin/time -v ~/projects/projects2014-metagenome/metagraph/build_test/metagraph transform_anno -v \
                --anno-type brwt \
                --greedy --relax-arity 20 \
                --parallel 30 \
                -o ~/metagenome/data/mantis/subsets/clean_rnaseq/mtg/brwt_${N}k \
                2>&1 | tee ~/metagenome/data/mantis/subsets/clean_rnaseq/mtg/logs/clean_primary_subset_${N}k_brwt.log"; \
    bsub -J "clean_brwt_${N}_relax" \
         -w "clean_brwt_${N}" \
         -oo ~/metagenome/data/mantis/subsets/clean_rnaseq/mtg/logs/clean_primary_subset_${N}k_brwt_relax.lsf \
         -W 48:00 \
         -n 15 -R "rusage[mem=41000] span[hosts=1]" \
        "/usr/bin/time -v ~/projects/projects2014-metagenome/metagraph/build_test/metagraph relax_brwt -v \
                --relax-arity 20 \
                --parallel 30 \
                -o ~/metagenome/data/mantis/subsets/clean_rnaseq/mtg/brwt_relaxed_${N}k \
                ~/metagenome/data/mantis/subsets/clean_rnaseq/mtg/brwt_${N}k.brwt.annodbg \
                2>&1 | tee ~/metagenome/data/mantis/subsets/clean_rnaseq/mtg/logs/clean_primary_subset_${N}k_brwt_relax.log"; \
done
```


## RefSeq (Fungi)

```bash

for N in {1..8}; do
    bsub -J "fungi_subset_${N}" \
         -oo ~/metagenome/data/mantis/subsets/refseq/mtg/logs/fungi_subset_${N}.lsf \
         -W 8:00 \
         -n 15 -R "rusage[mem=22000] span[hosts=1]" \
        "cat ~/metagenome/data/mantis/subsets/refseq/${N}.txt \
            | /usr/bin/time -v ~/projects/projects2014-metagenome/metagraph/build_test/metagraph build -v \
                -k 31 \
                --parallel 30 \
                --mem-cap-gb 300 \
                --disk-cap-gb 10000 \
                --disk-swap ~/metagenome/data/mantis/subsets/refseq/mtg/ \
                -o ~/metagenome/data/mantis/subsets/refseq/mtg/fungi_subset_${N} \
                2>&1"; \
done


for N in {1..8}; do
    mkdir ~/metagenome/data/mantis/subsets/refseq/mtg/cols_${N};
    bsub -J "fungi_annotate_${N}" \
         -oo ~/metagenome/data/mantis/subsets/refseq/mtg/logs/fungi_subset_${N}_annotate.lsf \
         -w "fungi_subset_${N}" \
         -W 15:00 \
         -n 30 -R "rusage[mem=15000] span[hosts=1]" \
        "cat ~/metagenome/data/mantis/subsets/refseq/${N}.txt \
            | /usr/bin/time -v ~/projects/projects2014-metagenome/metagraph/build_test/metagraph annotate -v \
                --parallel 60 \
                --anno-filename \
                --separately \
                -o ~/metagenome/data/mantis/subsets/refseq/mtg/cols_${N} \
                -i ~/metagenome/data/mantis/subsets/refseq/mtg/fungi_subset_${N}.dbg \
                2>&1 | tee ~/metagenome/data/mantis/subsets/refseq/mtg/logs/fungi_subset_${N}_annotate.log"; \
done


DIR=~/metagenome/data/row_diff/subsets/refseq/nobackup/mtg
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_master/metagraph
for DEPTH in {100,75,50,25,10}; do
    THREADS=72
    for N in {8,}; do
        NAME=refseq_rd_${DEPTH}d_${N}_${THREADS}t;
        RD_DIR=${DIR}/${NAME};
        rm -r ${RD_DIR};
        mkdir ${RD_DIR};
        ln -s ${DIR}/fungi_subset_${N}.dbg \
              ${RD_DIR}/graph.dbg;
        bsub -J "${NAME}" \
             -oo ${DIR}/logs/${NAME}.lsf \
             -W 10:00 \
             -n 36 -R "rusage[mem=9500] span[hosts=1]" \
            "find ${DIR}/cols_${N} -name \"*.column.annodbg\" \
                | /usr/bin/time -v $METAGRAPH transform_anno -v \
                    --anno-type row_diff \
                    --max-path-length $DEPTH \
                    --mem-cap-gb 300 \
                    --parallel $THREADS \
                    --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
                    -o ${RD_DIR}/out \
                    -i ${RD_DIR}/graph.dbg \
                    2>&1 | tee ${DIR}/logs/${NAME}.log";
        bsub -J "${NAME}_opt" \
             -w "${NAME}" \
             -oo ${DIR}/logs/${NAME}_opt.lsf \
             -W 10:00 \
             -n 36 -R "rusage[mem=9500] span[hosts=1]" \
            "find ${DIR}/cols_${N} -name \"*.column.annodbg\" \
                | /usr/bin/time -v $METAGRAPH transform_anno -v \
                    --anno-type row_diff \
                    --max-path-length $DEPTH \
                    --mem-cap-gb 600 \
                    --parallel $THREADS \
                    --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
                    --optimize \
                    -o ${RD_DIR}/out \
                    -i ${RD_DIR}/graph.dbg \
                    2>&1 | tee ${DIR}/logs/${NAME}_opt.log"; \
        bsub -J "${NAME}_multi_brwt" \
             -w "${NAME}_opt" \
             -oo ${DIR}/logs/${NAME}_multi_brwt.lsf \
             -W 10:00 \
             -n 36 -R "rusage[mem=9500] span[hosts=1]" \
            "find ${RD_DIR} -name \"*.row_diff.annodbg\" \
                | /usr/bin/time -v $METAGRAPH transform_anno -v \
                    --anno-type row_diff_brwt \
                    --greedy \
                    --parallel $THREADS \
                    --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
                    -o ${DIR}/multi_brwt/${NAME} \
                    -i ${RD_DIR}/graph.dbg \
                    2>&1 | tee ${DIR}/logs/${NAME}_multi_brwt.log"; \
        bsub -J "${NAME}_multi_brwt_relax" \
             -w "${NAME}_multi_brwt" \
             -oo ${DIR}/logs/${NAME}_multi_brwt_relax.lsf \
             -W 10:00 \
             -n 36 -R "rusage[mem=9500] span[hosts=1]" \
            "/usr/bin/time -v $METAGRAPH relax_brwt -v \
                    --relax-arity 20 \
                    --parallel $THREADS \
                    --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
                    -o ${DIR}/multi_brwt/${NAME}_relaxed \
                    ${DIR}/multi_brwt/${NAME}.row_diff_brwt.annodbg \
                    2>&1 | tee ${DIR}/logs/${NAME}_multi_brwt_relax.log"; \
        bsub -J "${NAME}_row_sparse" \
             -w "${NAME}_opt" \
             -oo ${DIR}/logs/${NAME}_row_sparse.lsf \
             -W 10:00 \
             -n 36 -R "rusage[mem=9500] span[hosts=1]" \
            "find ${RD_DIR} -name \"*.row_diff.annodbg\" \
                | /usr/bin/time -v $METAGRAPH transform_anno -v \
                    --anno-type row_diff_sparse \
                    --parallel $THREADS \
                    --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
                    -o ${DIR}/row_sparse/${NAME} \
                    -i ${RD_DIR}/graph.dbg \
                    2>&1 | tee ${DIR}/logs/${NAME}_row_sparse.log"; \
    done
done


for N in {1..8}; do
    bsub -J "fungi_rb_brwt_${N}" \
         -w "fungi_annotate_${N}" \
         -oo ~/metagenome/data/mantis/subsets/refseq/mtg/logs/fungi_subset_${N}_rb_brwt.lsf \
         -W 48:00 \
         -n 15 -R "rusage[mem=41000] span[hosts=1]" \
        "find ~/metagenome/data/mantis/subsets/refseq/mtg/cols_${N} -name \"*.column.annodbg\" \
            | /usr/bin/time -v ~/projects/projects2014-metagenome/metagraph/build_test/metagraph transform_anno -v \
                --anno-type rb_brwt \
                --greedy --relax-arity 20 \
                --parallel 30 \
                -o ~/metagenome/data/mantis/subsets/refseq/mtg/rb_brwt_${N} \
                2>&1 | tee ~/metagenome/data/mantis/subsets/refseq/mtg/logs/fungi_subset_${N}_rb_brwt.log"; \
done


for N in {1..8}; do
    bsub -J "fungi_brwt_${N}" \
         -w "fungi_annotate_${N}" \
         -oo ~/metagenome/data/mantis/subsets/refseq/mtg/logs/fungi_subset_${N}_brwt.lsf \
         -W 48:00 \
         -n 15 -R "rusage[mem=11000] span[hosts=1]" \
        "find ~/metagenome/data/mantis/subsets/refseq/mtg/cols_${N} -name \"*.column.annodbg\" \
            | /usr/bin/time -v ~/projects/projects2014-metagenome/metagraph/build_test/metagraph transform_anno -v \
                --anno-type brwt \
                --greedy --relax-arity 20 \
                --parallel 30 \
                -o ~/metagenome/data/mantis/subsets/refseq/mtg/brwt_${N} \
                2>&1 | tee ~/metagenome/data/mantis/subsets/refseq/mtg/logs/fungi_subset_${N}_brwt.log"; \
    bsub -J "fungi_brwt_${N}_relax" \
         -w "fungi_brwt_${N}" \
         -oo ~/metagenome/data/mantis/subsets/refseq/mtg/logs/fungi_subset_${N}_brwt_relax.lsf \
         -W 10:00 \
         -n 15 -R "rusage[mem=11000] span[hosts=1]" \
        "/usr/bin/time -v ~/projects/projects2014-metagenome/metagraph/build_test/metagraph relax_brwt -v \
                --relax-arity 20 \
                --parallel 30 \
                -o ~/metagenome/data/mantis/subsets/refseq/mtg/brwt_relaxed_${N} \
                ~/metagenome/data/mantis/subsets/refseq/mtg/brwt_${N}.brwt.annodbg \
                2>&1 | tee ~/metagenome/data/mantis/subsets/refseq/mtg/logs/fungi_subset_${N}_brwt_relax.log"; \
done
```


## Query scripts

```bash
DIR=~/metagenome/data/row_diff/subsets/rnaseq/nobackup/mtg
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_test/metagraph

for QUERY in {transcripts_100.fa,transcripts_1000.fa}; do
    for DEPTH in {10,25,50,75,100}; do
        ANNO=rd_${DEPTH}d_10k_72t.row_diff_sparse.annodbg
        bsub -J "query_${ANNO}_${QUERY}" \
             -oo ~/metagenome/data/row_diff/query/${ANNO}_${QUERY}.lsf \
             -W 4:00 \
             -n 36 -R "rusage[mem=9500] span[hosts=1] select[model==XeonGold_6140]" \
                  "/usr/bin/time -v $METAGRAPH query -v \
                          --count-labels --fast \
                          -i $DIR/graph_primary_subset_10k_primary.dbg \
                          -a $DIR/row_sparse/$ANNO \
                          ~/metagenome/data/row_diff/query/$QUERY \
                          > ~/metagenome/data/row_diff/query/${ANNO}_${QUERY}.out \
                          2> ~/metagenome/data/row_diff/query/${ANNO}_${QUERY}.log"
    done
done

for QUERY in {transcripts_100.fa,transcripts_1000.fa}; do
    for DEPTH in {10,25,50,75,100}; do
        ANNO=rd_${DEPTH}d_10k_72t_relaxed.row_diff_brwt.annodbg
        bsub -J "query_${ANNO}_${QUERY}" \
         -oo ~/metagenome/data/row_diff/query/${ANNO}_${QUERY}.lsf \
         -W 4:00 \
         -n 36 -R "rusage[mem=9500] span[hosts=1] select[model==XeonGold_6140]" \
            "/usr/bin/time -v $METAGRAPH query -v \
                --count-labels --fast \
                -i $DIR/graph_primary_subset_10k_primary.dbg \
                -a $DIR/multi_brwt/$ANNO \
                ~/metagenome/data/row_diff/query/$QUERY \
                > ~/metagenome/data/row_diff/query/${ANNO}_${QUERY}.out \
                2> ~/metagenome/data/row_diff/query/${ANNO}_${QUERY}.log"
    done
done

for QUERY in {transcripts_100.fa,transcripts_1000.fa}; do
    ANNO=brwt_relaxed_10k.brwt.annodbg
    bsub -J "query_${ANNO}_${QUERY}" \
         -oo ~/metagenome/data/row_diff/query/${ANNO}_${QUERY}.lsf \
         -W 4:00 \
         -n 36 -R "rusage[mem=9500] span[hosts=1] select[model==XeonGold_6140]" \
            "/usr/bin/time -v $METAGRAPH query -v \
                    --count-labels --fast \
                    -i $DIR/graph_primary_subset_10k_primary.dbg \
                    -a $DIR/multi_brwt/$ANNO \
                    ~/metagenome/data/row_diff/query/$QUERY \
                    > ~/metagenome/data/row_diff/query/${ANNO}_${QUERY}.out \
                    2> ~/metagenome/data/row_diff/query/${ANNO}_${QUERY}.log"
done
```
