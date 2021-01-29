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


DEPTH=100
THREADS=72
DIR=~/metagenome/data/mantis/subsets/mtg
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_release/metagraph
for N in {1..10}; do
    NAME=rd_${DEPTH}d_${N}k_${THREADS}t;
    RD_DIR=${DIR}/${NAME};
    rm -r ${RD_DIR};
    mkdir ${RD_DIR};
    ln -s ${DIR}/graph_primary_subset_${N}k.dbg \
          ${RD_DIR}/graph.dbg;
    bsub -J "${NAME}" \
         -oo ${DIR}/logs/${NAME}.lsf \
         -W 10:00 \
         -n 36 -R "rusage[mem=19000] span[hosts=1] select[model==XeonGold_6140]" \
        "find ${DIR}/cols_${N}k -name \"*.column.annodbg\" \
            | /usr/bin/time -v $METAGRAPH transform_anno -v \
                --anno-type row_diff \
                --max-path-length $DEPTH \
                --mem-cap-gb 600 \
                --parallel $THREADS \
                -o ${RD_DIR}/out \
                -i ${RD_DIR}/graph.dbg \
                2>&1 | tee ${DIR}/logs/${NAME}.log";
    bsub -J "${NAME}_opt" \
         -w "${NAME}" \
         -oo ${DIR}/logs/${NAME}_opt.lsf \
         -W 10:00 \
         -n 36 -R "rusage[mem=19000] span[hosts=1] select[model==XeonGold_6140]" \
        "find ${DIR}/cols_${N}k -name \"*.column.annodbg\" \
            | /usr/bin/time -v $METAGRAPH transform_anno -v \
                --anno-type row_diff \
                --max-path-length $DEPTH \
                --mem-cap-gb 600 \
                --parallel $THREADS \
                --optimize \
                -o ${RD_DIR}/out \
                -i ${RD_DIR}/graph.dbg \
                2>&1 | tee ${DIR}/logs/${NAME}_opt.log"; \
    bsub -J "${NAME}_multi_brwt" \
         -w "${NAME}_opt" \
         -oo ${DIR}/logs/${NAME}_multi_brwt.lsf \
         -W 10:00 \
         -n 36 -R "rusage[mem=19000] span[hosts=1] select[model==XeonGold_6140]" \
        "find ${RD_DIR} -name \"*.row_diff.annodbg\" \
            | /usr/bin/time -v $METAGRAPH transform_anno -v \
                --anno-type row_diff_brwt \
                --greedy \
                --parallel $THREADS \
                -o ${DIR}/multi_brwt/${NAME} \
                --anchors-file ${RD_DIR}/graph.dbg.anchors \
                2>&1 | tee ${DIR}/logs/${NAME}_multi_brwt.log"; \
    bsub -J "${NAME}_multi_brwt_relax" \
         -w "${NAME}_multi_brwt" \
         -oo ${DIR}/logs/${NAME}_multi_brwt_relax.lsf \
         -W 10:00 \
         -n 36 -R "rusage[mem=19000] span[hosts=1] select[model==XeonGold_6140]" \
        "/usr/bin/time -v $METAGRAPH relax_brwt -v \
                --relax-arity 20 \
                --parallel $THREADS \
                -o ${DIR}/multi_brwt/${NAME}_relaxed \
                ${DIR}/multi_brwt/${NAME}.row_diff_brwt.annodbg \
                2>&1 | tee ${DIR}/logs/${NAME}_multi_brwt_relax.log"; \
done


DEPTH=100
for N in {1..10}; do
    bsub -J "rd_${DEPTH}_${N}_to_row_sparse" \
         -oo ~/metagenome/data/mantis/subsets/mtg/logs/graph_primary_subset_${N}k_rd_${DEPTH}_to_row_sparse.lsf \
         -W 20:00 \
         -n 15 -R "rusage[mem=51000] span[hosts=1]" \
        "find ~/metagenome/data/mantis/subsets/mtg/rd_${DEPTH}_${N}k -name \"*.row_diff.annodbg\" \
            | /usr/bin/time -v ~/projects/projects2014-metagenome/metagraph/build_test/metagraph transform_anno -v \
                --anno-type row_diff_sparse \
                --parallel 30 \
                -o ~/metagenome/data/mantis/subsets/mtg/row_sparse/rd_${DEPTH}_${N}k \
                --anchors-file ~/metagenome/data/mantis/subsets/mtg/rd_${DEPTH}_${N}k/graph.dbg.anchors \
                2>&1 | tee ~/metagenome/data/mantis/subsets/mtg/logs/graph_primary_subset_${N}k_rd_${DEPTH}_to_row_sparse.log"; \
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
done


for N in {1..10}; do
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


DEPTH=100
for N in {1..10}; do
    rm -r ~/metagenome/data/mantis/subsets/clean_rnaseq/mtg/rd_${DEPTH}_${N}k;
    mkdir ~/metagenome/data/mantis/subsets/clean_rnaseq/mtg/rd_${DEPTH}_${N}k;
    ln -s ~/metagenome/data/mantis/subsets/clean_rnaseq/mtg/clean_primary_subset_${N}k.dbg \
          ~/metagenome/data/mantis/subsets/clean_rnaseq/mtg/rd_${DEPTH}_${N}k/graph.dbg;
    bsub -J "clean_rd_${DEPTH}_${N}" \
         -oo ~/metagenome/data/mantis/subsets/clean_rnaseq/mtg/logs/clean_primary_subset_${N}k_rd_${DEPTH}.lsf \
         -W 20:00 \
         -n 15 -R "rusage[mem=40000] span[hosts=1]" \
        "find ~/metagenome/data/mantis/subsets/clean_rnaseq/mtg/cols_${N}k -name \"*.column.annodbg\" \
            | /usr/bin/time -v ~/projects/projects2014-metagenome/metagraph/build_test/metagraph transform_anno -v \
                --anno-type row_diff \
                --max-path-length $DEPTH \
                --mem-cap-gb 500 \
                --parallel 30 \
                -o ~/metagenome/data/mantis/subsets/clean_rnaseq/mtg/rd_${DEPTH}_${N}k/out \
                -i ~/metagenome/data/mantis/subsets/clean_rnaseq/mtg/rd_${DEPTH}_${N}k/graph.dbg \
                2>&1 | tee ~/metagenome/data/mantis/subsets/clean_rnaseq/mtg/logs/clean_primary_subset_${N}k_rd_${DEPTH}.log"; \
done


DEPTH=100
for N in {1..10}; do
    bsub -J "clean_rd_${DEPTH}_${N}_opt" \
         -w "clean_rd_${DEPTH}_${N}" \
         -oo ~/metagenome/data/mantis/subsets/clean_rnaseq/mtg/logs/clean_primary_subset_${N}k_rd_${DEPTH}_opt.lsf \
         -W 20:00 \
         -n 15 -R "rusage[mem=40000] span[hosts=1]" \
        "find ~/metagenome/data/mantis/subsets/clean_rnaseq/mtg/cols_${N}k -name \"*.column.annodbg\" \
            | /usr/bin/time -v ~/projects/projects2014-metagenome/metagraph/build_test/metagraph transform_anno -v \
                --anno-type row_diff \
                --max-path-length $DEPTH \
                --mem-cap-gb 500 \
                --parallel 30 \
                --optimize \
                -o ~/metagenome/data/mantis/subsets/clean_rnaseq/mtg/rd_${DEPTH}_${N}k/out \
                -i ~/metagenome/data/mantis/subsets/clean_rnaseq/mtg/rd_${DEPTH}_${N}k/graph.dbg \
                2>&1 | tee ~/metagenome/data/mantis/subsets/clean_rnaseq/mtg/logs/clean_primary_subset_${N}k_rd_${DEPTH}_opt.log"; \
done


mkdir ~/metagenome/data/mantis/subsets/clean_rnaseq/mtg/row_sparse
DEPTH=100
for N in {1..10}; do
    bsub -J "clean_rd_${DEPTH}_${N}_row_sparse" \
         -w "clean_rd_${DEPTH}_${N}_opt" \
         -oo ~/metagenome/data/mantis/subsets/clean_rnaseq/mtg/logs/clean_primary_subset_${N}k_rd_${DEPTH}_row_sparse.lsf \
         -W 20:00 \
         -n 15 -R "rusage[mem=41000] span[hosts=1]" \
        "find ~/metagenome/data/mantis/subsets/clean_rnaseq/mtg/rd_${DEPTH}_${N}k -name \"*.row_diff.annodbg\" \
            | /usr/bin/time -v ~/projects/projects2014-metagenome/metagraph/build_test/metagraph transform_anno -v \
                --anno-type row_diff_sparse \
                --parallel 30 \
                -o ~/metagenome/data/mantis/subsets/clean_rnaseq/mtg/row_sparse/rd_${DEPTH}_${N}k \
                --anchors-file ~/metagenome/data/mantis/subsets/clean_rnaseq/mtg/rd_${DEPTH}_${N}k/graph.dbg.anchors \
                2>&1 | tee ~/metagenome/data/mantis/subsets/clean_rnaseq/mtg/logs/clean_primary_subset_${N}k_rd_${DEPTH}_row_sparse.log"; \
done


mkdir ~/metagenome/data/mantis/subsets/clean_rnaseq/mtg/multi_brwt
DEPTH=100
for N in {1..10}; do
    bsub -J "clean_rd_${DEPTH}_${N}_multi_brwt" \
         -w "clean_rd_${DEPTH}_${N}_opt" \
         -oo ~/metagenome/data/mantis/subsets/clean_rnaseq/mtg/logs/clean_primary_subset_${N}k_rd_${DEPTH}_multi_brwt.lsf \
         -W 20:00 \
         -n 15 -R "rusage[mem=41000] span[hosts=1]" \
        "find ~/metagenome/data/mantis/subsets/clean_rnaseq/mtg/rd_${DEPTH}_${N}k -name \"*.row_diff.annodbg\" \
            | /usr/bin/time -v ~/projects/projects2014-metagenome/metagraph/build_test/metagraph transform_anno -v \
                --anno-type row_diff_brwt \
                --greedy --relax-arity 20 \
                --parallel 30 \
                -o ~/metagenome/data/mantis/subsets/clean_rnaseq/mtg/multi_brwt/rd_${DEPTH}_${N}k \
                --anchors-file ~/metagenome/data/mantis/subsets/clean_rnaseq/mtg/rd_${DEPTH}_${N}k/graph.dbg.anchors \
                2>&1 | tee ~/metagenome/data/mantis/subsets/clean_rnaseq/mtg/logs/clean_primary_subset_${N}k_rd_${DEPTH}_multi_brwt.log"; \
done


DEPTH=100
for N in {1..10}; do
    bsub -J "clean_rd_${DEPTH}_${N}_multi_brwt_relax" \
         -w "clean_rd_${DEPTH}_${N}_multi_brwt" \
         -oo ~/metagenome/data/mantis/subsets/clean_rnaseq/mtg/logs/clean_primary_subset_${N}k_rd_${DEPTH}_multi_brwt_relax.lsf \
         -W 20:00 \
         -n 15 -R "rusage[mem=11000] span[hosts=1]" \
        "/usr/bin/time -v ~/projects/projects2014-metagenome/metagraph/build_test/metagraph relax_brwt -v \
                --relax-arity 20 \
                --parallel 30 \
                -o ~/metagenome/data/mantis/subsets/clean_rnaseq/mtg/multi_brwt/rd_${DEPTH}_${N}k_relaxed \
                ~/metagenome/data/mantis/subsets/clean_rnaseq/mtg/multi_brwt/rd_${DEPTH}_${N}k.row_diff_brwt.annodbg \
                2>&1 | tee ~/metagenome/data/mantis/subsets/clean_rnaseq/mtg/logs/clean_primary_subset_${N}k_rd_${DEPTH}_multi_brwt_relax.log"; \
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
done


for N in {1..10}; do
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


DEPTH=100
for N in {1..8}; do
    mkdir ~/metagenome/data/mantis/subsets/refseq/mtg/rd_${DEPTH}_${N};
    ln -s ~/metagenome/data/mantis/subsets/refseq/mtg/fungi_subset_${N}.dbg \
          ~/metagenome/data/mantis/subsets/refseq/mtg/rd_${DEPTH}_${N}/graph.dbg;
    bsub -J "fungi_rd_${DEPTH}_${N}" \
         -w "fungi_annotate_${N}" \
         -oo ~/metagenome/data/mantis/subsets/refseq/mtg/logs/fungi_subset_${N}_rd_${DEPTH}.lsf \
         -W 20:00 \
         -n 15 -R "rusage[mem=21000] span[hosts=1]" \
        "find ~/metagenome/data/mantis/subsets/refseq/mtg/cols_${N} -name \"*.column.annodbg\" \
            | /usr/bin/time -v ~/projects/projects2014-metagenome/metagraph/build_test/metagraph transform_anno -v \
                --anno-type row_diff \
                --max-path-length $DEPTH \
                --mem-cap-gb 300 \
                --parallel 30 \
                -o ~/metagenome/data/mantis/subsets/refseq/mtg/rd_${DEPTH}_${N}/out \
                -i ~/metagenome/data/mantis/subsets/refseq/mtg/rd_${DEPTH}_${N}/graph.dbg \
                2>&1 | tee ~/metagenome/data/mantis/subsets/refseq/mtg/logs/fungi_subset_${N}_rd_${DEPTH}.log"; \
done


DEPTH=100
for N in {1..8}; do
    bsub -J "fungi_rd_${DEPTH}_${N}_opt" \
         -w "fungi_rd_${DEPTH}_${N}" \
         -oo ~/metagenome/data/mantis/subsets/refseq/mtg/logs/fungi_subset_${N}_rd_${DEPTH}_opt.lsf \
         -W 20:00 \
         -n 15 -R "rusage[mem=21000] span[hosts=1]" \
        "find ~/metagenome/data/mantis/subsets/refseq/mtg/cols_${N} -name \"*.column.annodbg\" \
            | /usr/bin/time -v ~/projects/projects2014-metagenome/metagraph/build_test/metagraph transform_anno -v \
                --anno-type row_diff \
                --max-path-length $DEPTH \
                --mem-cap-gb 300 \
                --parallel 30 \
                --optimize \
                -o ~/metagenome/data/mantis/subsets/refseq/mtg/rd_${DEPTH}_${N}/out \
                -i ~/metagenome/data/mantis/subsets/refseq/mtg/rd_${DEPTH}_${N}/graph.dbg \
                2>&1 | tee ~/metagenome/data/mantis/subsets/refseq/mtg/logs/fungi_subset_${N}_rd_${DEPTH}_opt.log"; \
done


mkdir ~/metagenome/data/mantis/subsets/refseq/mtg/row_sparse
DEPTH=100
for N in {1..8}; do
    bsub -J "fungi_rd_${DEPTH}_${N}_row_sparse" \
         -w "fungi_rd_${DEPTH}_${N}_opt" \
         -oo ~/metagenome/data/mantis/subsets/refseq/mtg/logs/fungi_subset_${N}_rd_${DEPTH}_row_sparse.lsf \
         -W 20:00 \
         -n 15 -R "rusage[mem=41000] span[hosts=1]" \
        "find ~/metagenome/data/mantis/subsets/refseq/mtg/rd_${DEPTH}_${N} -name \"*.row_diff.annodbg\" \
            | /usr/bin/time -v ~/projects/projects2014-metagenome/metagraph/build_test/metagraph transform_anno -v \
                --anno-type row_diff_sparse \
                --parallel 30 \
                -o ~/metagenome/data/mantis/subsets/refseq/mtg/row_sparse/rd_${DEPTH}_${N} \
                --anchors-file ~/metagenome/data/mantis/subsets/refseq/mtg/rd_${DEPTH}_${N}/graph.dbg.anchors \
                2>&1 | tee ~/metagenome/data/mantis/subsets/refseq/mtg/logs/fungi_subset_${N}_rd_${DEPTH}_row_sparse.log"; \
done


mkdir ~/metagenome/data/mantis/subsets/refseq/mtg/multi_brwt
DEPTH=100
for N in {1..8}; do
    bsub -J "fungi_rd_${DEPTH}_${N}_multi_brwt" \
         -w "fungi_rd_${DEPTH}_${N}_opt" \
         -oo ~/metagenome/data/mantis/subsets/refseq/mtg/logs/fungi_subset_${N}_rd_${DEPTH}_multi_brwt.lsf \
         -W 20:00 \
         -n 15 -R "rusage[mem=41000] span[hosts=1]" \
        "find ~/metagenome/data/mantis/subsets/refseq/mtg/rd_${DEPTH}_${N} -name \"*.row_diff.annodbg\" \
            | /usr/bin/time -v ~/projects/projects2014-metagenome/metagraph/build_test/metagraph transform_anno -v \
                --anno-type row_diff_brwt \
                --greedy --relax-arity 20 \
                --parallel 30 \
                -o ~/metagenome/data/mantis/subsets/refseq/mtg/multi_brwt/rd_${DEPTH}_${N} \
                --anchors-file ~/metagenome/data/mantis/subsets/refseq/mtg/rd_${DEPTH}_${N}/graph.dbg.anchors \
                2>&1 | tee ~/metagenome/data/mantis/subsets/refseq/mtg/logs/fungi_subset_${N}_rd_${DEPTH}_multi_brwt.log"; \
done


DEPTH=100
for N in {1..8}; do
    bsub -J "fungi_rd_${DEPTH}_${N}_multi_brwt_relax" \
         -w "fungi_rd_${DEPTH}_${N}_multi_brwt" \
         -oo ~/metagenome/data/mantis/subsets/refseq/mtg/logs/fungi_subset_${N}_rd_${DEPTH}_multi_brwt_relax.lsf \
         -W 20:00 \
         -n 15 -R "rusage[mem=11000] span[hosts=1]" \
        "/usr/bin/time -v ~/projects/projects2014-metagenome/metagraph/build_test/metagraph relax_brwt -v \
                --relax-arity 20 \
                --parallel 30 \
                -o ~/metagenome/data/mantis/subsets/refseq/mtg/multi_brwt/rd_${DEPTH}_${N}_relaxed \
                ~/metagenome/data/mantis/subsets/refseq/mtg/multi_brwt/rd_${DEPTH}_${N}.row_diff_brwt.annodbg \
                2>&1 | tee ~/metagenome/data/mantis/subsets/refseq/mtg/logs/fungi_subset_${N}_rd_${DEPTH}_multi_brwt_relax.log"; \
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
done


for N in {1..8}; do
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
