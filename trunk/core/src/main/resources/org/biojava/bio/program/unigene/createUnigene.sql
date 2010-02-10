# unigene schema
#

# The unigene entry
CREATE TABLE cluster (
  id varchar(10) unique NOT NULL,
  title varchar(255),
  gene varchar(255),
  cytoband varchar(7),
  locuslink int,
  gnm_terminus char(1),
  chromosome char(2),
  unique key cluster_id (id)
);

# Expression list for a unigene entry
CREATE TABLE express (
  cluster_id varchar(10) not null,
  tissue varchar(255),
  unique key express_key (cluster_id, tissue)
);

# sts for a unigene entry
CREATE TABLE sts (
  cluster_id varchar(10) not null,
  acc varchar(10),
  name varchar(10),
  unist varchar(10),
  key sts_key (cluster_id)
);

# TXMAP for a unigene entry
CREATE TABLE txmap (
  cluster_id varchar(10) not null,
  id varchar(16),
  key txmap_key (cluster_id)
);

# protsim for a unigene entry
CREATE TABLE protsim (
  cluster_id varchar(10) not null,
  org varchar(255),
  protgi varchar(10),
  protid varchar(15),
  pct int(4),
  aln int(4),
  key protsim_key (cluster_id)
);

# sequences in a unigene cluster
CREATE TABLE sequence (
  cluster_id varchar(10) not null,
  acc varchar(15),
  nid varchar(15),
  clone varchar(15),
  clone_end char(2),
  lid int(5),
  pid varchar(15),
  key sequence_key (cluster_id)
);

