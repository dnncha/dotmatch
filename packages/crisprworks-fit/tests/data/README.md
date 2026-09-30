# Upstream regression fixture

`mageck2-counts.tsv` is copied without changes from MAGeCK2's public smoke-test
count table at commit `630aea0b6fc152a81006435f21d629297275c911`:

https://github.com/davidliwei/mageck2/blob/630aea0b6fc152a81006435f21d629297275c911/tests/data/count_table.txt

Git blob SHA: `633f72aaeb99a0b2e7bd40e4742056f0d94340ed`.

It contains 999 guides across 100 genes and four samples, including zero
counts. `design.tsv` assigns the two initial samples to the shared baseline and
the two final samples to separate effects, as in upstream's MLE regression
tests. MAGeCK2's BSD-3-Clause notice is retained in the package LICENSE.

This checks compatibility with the upstream fixture. It does not establish
biological accuracy or performance across genome-wide industry datasets.
