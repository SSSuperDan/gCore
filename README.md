# gCore: Exploring Cross-layer Cohesiveness in Multi-layer graphs


## Execution

### Usage

```commandline
gCore mode [graph_path graph_name] [options]
```
#### Mode:

   	cs      Run case study.
	grk     Generate random \mathbf{k} vectors.
	grp     Generate random \mathbf{p} vectors.
	grkp    Generate random (\mathbf{k},\mathbf{p}) vector pairs.
	pgcs	Compare the efficiency of CORE, dCC and GCS on pillar multi-layer graphs.
	bp      Build P-tree.
    bkp     Build KP-tree.
	gcs     Compare the efficiency of CORE, dCC(pillar mlg only), RCD(general mlg only), GCS and (K)P-tree-based search (under different P-tree compaction levels).
	gcii	Show statistics of GCI, including construction time, memory cost and number of nodes.
	kv  	Compare the k-value of nodes in the (\mathbf{k},\mathbf{p})-core, the (\mathbf{k},\mathbf{ck})-rcd and the k-core.
	pv  	Compare the p-value of nodes in the (\mathbf{k},\mathbf{p})-core, the (\mathbf{k},\mathbf{ck})-rcd and the k-core.
	sm      Compute the size distribution when varying \mathbf{k}[i] and \mathbf{p}[i] for each layer i.
	smk     Compute the size distribution when fixing \mathbf{k} and varying \mathbf{p}[i] for each layer i.
	info	Compute basic information of the loaded multi-layer graph.
	
#### Options:
    -ntc	Number of testcases.
	-skf	File of sampled coreness vectors (\mathbf{k}).
	-spf	File of sampled fraction vectors (\mathbf{p}).
	-skpf	File of sampled (\mathbf{k}, \mathbf{p}) vector pairs.
	-k      Coreness vector (\mathbf{k}).
	-p      Neighbor coverage fraction vector (\mathbf{p}).
	-b/pb   P-tree builder.
	-s      Incremental step of coreness vectors to construct P-trees.
	-o      Output path.
	-ik     Start/Initial coreness vector.
	-ck     Cross-layer degree threshold vector.
	-pk     Coreness threshold for the selected layer.
	-ptf	P-tree file.
	-kptf	KP-tree file.
	-f2if	Fraction to index map file.
	-ps     Incremental step of fraction vectors to compute size matrices.
	-ek     End coreness vector.
	-d/dim  Dimension.
	-g      Graph to perform case study, "dblp" and "twitter" are availiable.

#### P-tree builder options:
    naive	Build P-tree without optimization.
    ne      Build compact P+-tree with subtree elimination.
    se      Build compact P+-DAG with subtree merge
    nese	Build compact P+-DAG with both subtree elimination and subtree merge.

