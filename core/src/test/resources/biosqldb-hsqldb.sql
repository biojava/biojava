-- $Id$
--
-- biosqldb-hsqldb.sql

-- Authors: Ewan Birney, Elia Stupka
-- Contributors: Hilmar Lapp, Aaron Mackey
--
-- Copyright Ewan Birney. You may use, modify, and distribute this code under
-- the same terms as Perl. See the Perl Artistic License.
--
-- comments to biosql - biosql-l@open-bio.org
--
-- Migration of the MySQL schema to InnoDB by Hilmar Lapp <hlapp at gmx.net>
-- Post-Cape Town changes by Hilmar Lapp.
-- Singapore changes by Hilmar Lapp and Aaron Mackey.
-- Migration to HSQLDB by Len Trigg <len at reeltwo.com>


-- See biosql schema documentation for general documentation regarding the 
-- schema. This file contains documetation specific to the hsqldb schema.

-- HSQLDB Version compatibility notes:
-- HSQLDB 1.7.1 has problems with null values in columns with UNIQUE constraints, several of these
-- constraints have been commented out for compatibility.
-- HSQLDB 1.7.2 alpha N. has problem with PreparedStatements that affect the BioJava binding (these 
-- will apparently be addressed by the 1.7.2 release)

CREATE TABLE biodatabase (
  	biodatabase_id 	INT NOT NULL IDENTITY,
  	name           	VARCHAR(128) NOT NULL,
	authority	VARCHAR(128),
	description	LONGVARCHAR,
        PRIMARY KEY (biodatabase_id),
  	UNIQUE (name)
);

CREATE INDEX db_auth on biodatabase(authority);


CREATE TABLE taxon (
       taxon_id		INT NOT NULL IDENTITY,
       ncbi_taxon_id 	INT,
       parent_taxon_id	INT,
       node_rank	VARCHAR(32),
       genetic_code	TINYINT,
       mito_genetic_code TINYINT,
       left_value	INT,
       right_value	INT,
       PRIMARY KEY (taxon_id),
       UNIQUE (ncbi_taxon_id)
);
-- HSQLDB 1.7.1 UNIQUE BUG
--       UNIQUE (left_value),
--       UNIQUE (right_value)

CREATE INDEX taxparent ON taxon(parent_taxon_id);


CREATE TABLE taxon_name (
       taxon_id		INT NOT NULL,
       name		VARCHAR(255) NOT NULL,
       name_class	VARCHAR(32) NOT NULL,
       UNIQUE (taxon_id,name,name_class)
);

CREATE INDEX taxnametaxonid ON taxon_name(taxon_id);
CREATE INDEX taxnamename    ON taxon_name(name);


CREATE TABLE ontology (
       	ontology_id        INT NOT NULL IDENTITY,
       	name	   	   VARCHAR(32) NOT NULL,
       	definition	   LONGVARCHAR,
        PRIMARY KEY (ontology_id),
	UNIQUE (name)
);


CREATE TABLE term (
       	term_id   INT NOT NULL IDENTITY,
       	name	   	   VARCHAR(255) NOT NULL,
       	definition	   LONGVARCHAR,
	identifier	   VARCHAR(40),
	is_obsolete	   CHAR(1),
	ontology_id	   INT NOT NULL,
        PRIMARY KEY (term_id),
	UNIQUE (name,ontology_id)
);
-- HSQLDB 1.7.1 UNIQUE BUG
--	UNIQUE (identifier)

CREATE INDEX term_ont ON term(ontology_id);


-- We should use the field name "name" instead of "synonym" 
-- (because it is a reserved word in some RDBMS)
CREATE TABLE term_synonym (
       name		  VARCHAR(255) NOT NULL,
       term_id		  INT NOT NULL,
       PRIMARY KEY (term_id,name)
);


CREATE TABLE term_dbxref (
       	term_id	          INT NOT NULL,
       	dbxref_id         INT NOT NULL,
	rank		  SMALLINT,
	PRIMARY KEY (term_id, dbxref_id)
);

CREATE INDEX trmdbxref_dbxrefid ON term_dbxref(dbxref_id);


CREATE TABLE term_relationship (
        term_relationship_id INT NOT NULL IDENTITY,
       	subject_term_id	INT NOT NULL,
       	predicate_term_id    INT NOT NULL,
       	object_term_id       INT NOT NULL,
	ontology_id	INT NOT NULL,
        PRIMARY KEY (term_relationship_id),
	UNIQUE (subject_term_id,predicate_term_id,object_term_id,ontology_id)
);

CREATE INDEX trmrel_predicateid ON term_relationship(predicate_term_id);
CREATE INDEX trmrel_objectid ON term_relationship(object_term_id);
CREATE INDEX trmrel_ontid ON term_relationship(ontology_id);


-- This lets one associate a single term with a term_relationship 
-- effecively allowing us to treat triples as 1st class terms.
-- http://www.open-bio.org/pipermail/biosql-l/2003-October/000455.html
CREATE TABLE term_relationship_term (
        term_relationship_id INT NOT NULL,
        term_id              INT NOT NULL,
        UNIQUE ( term_relationship_id ),
        UNIQUE ( term_id ) 
);

CREATE TABLE term_path (
        term_path_id         INT NOT NULL IDENTITY,
       	subject_term_id	     INT NOT NULL,
       	predicate_term_id    INT NOT NULL,
       	object_term_id       INT NOT NULL,
	ontology_id          INT NOT NULL,
	distance	     INT,
        PRIMARY KEY (term_path_id),
	UNIQUE (subject_term_id,predicate_term_id,object_term_id,ontology_id,distance)
);

CREATE INDEX trmpath_predicateid ON term_path(predicate_term_id);
CREATE INDEX trmpath_objectid ON term_path(object_term_id);
CREATE INDEX trmpath_ontid ON term_path(ontology_id);


CREATE TABLE bioentry (
	bioentry_id	INT NOT NULL IDENTITY,
  	biodatabase_id  INT NOT NULL,
  	taxon_id     	INT,
  	name		VARCHAR(40) NOT NULL,
  	accession    	VARCHAR(40) NOT NULL,
  	identifier   	VARCHAR(40),
	division	VARCHAR(6),
  	description  	LONGVARCHAR,
  	version 	SMALLINT NOT NULL, 
        PRIMARY KEY (bioentry_id),
  	UNIQUE (accession,biodatabase_id,version),
        UNIQUE (identifier,biodatabase_id)
);
-- HSQLDB 1.7.1 UNIQUE BUG
--  	UNIQUE (identifier)

CREATE INDEX bioentry_name ON bioentry(name);
CREATE INDEX bioentry_db   ON bioentry(biodatabase_id);
CREATE INDEX bioentry_tax  ON bioentry(taxon_id);


CREATE TABLE bioentry_relationship (
        bioentry_relationship_id INT NOT NULL IDENTITY,
   	object_bioentry_id 	INT NOT NULL,
   	subject_bioentry_id 	INT NOT NULL,
   	term_id 		INT NOT NULL,
   	rank 			INT,
        PRIMARY KEY (bioentry_relationship_id),
	UNIQUE (object_bioentry_id,subject_bioentry_id,term_id)
);

CREATE INDEX bioentryrel_trm   ON bioentry_relationship(term_id);
CREATE INDEX bioentryrel_child ON bioentry_relationship(subject_bioentry_id);


CREATE TABLE bioentry_path (
   	object_bioentry_id 	INT NOT NULL,
   	subject_bioentry_id 	INT NOT NULL,
   	term_id 		INT NOT NULL,
	distance	     	INT,
	UNIQUE (object_bioentry_id,subject_bioentry_id,term_id,distance)
);

CREATE INDEX bioentrypath_trm   ON bioentry_path(term_id);
CREATE INDEX bioentrypath_child ON bioentry_path(subject_bioentry_id);


CREATE TABLE biosequence (
  	bioentry_id     INT NOT NULL,
  	version     	SMALLINT, 
  	length      	INT,
  	alphabet        VARCHAR(10),
  	seq 		LONGVARCHAR,
	PRIMARY KEY (bioentry_id)
);

-- add these only if you want them:
-- ALTER TABLE biosequence ADD COLUMN ( isoelec_pt NUMERIC(4,2) );
-- ALTER TABLE biosequence ADD COLUMN (	mol_wgt DOUBLE PRECISION );
-- ALTER TABLE biosequence ADD COLUMN ( perc_gc DOUBLE PRECISION );


CREATE TABLE dbxref (
        dbxref_id	INT NOT NULL IDENTITY,
        dbname          VARCHAR(40) NOT NULL,
        accession       VARCHAR(40) NOT NULL,
	version		SMALLINT NOT NULL,
        PRIMARY KEY (dbxref_id),
        UNIQUE(accession, dbname, version)
);

CREATE INDEX dbxref_db  ON dbxref(dbname);


CREATE TABLE dbxref_qualifier_value (
       	dbxref_id 		INT NOT NULL,
       	term_id 		INT NOT NULL,
  	rank  		   	SMALLINT DEFAULT 0 NOT NULL,
       	value			LONGVARCHAR,
	PRIMARY KEY (dbxref_id,term_id,rank)
);

CREATE INDEX dbxrefqual_dbx ON dbxref_qualifier_value(dbxref_id);
CREATE INDEX dbxrefqual_trm ON dbxref_qualifier_value(term_id);


CREATE TABLE bioentry_dbxref ( 
       	bioentry_id        INT NOT NULL,
       	dbxref_id          INT NOT NULL,
  	rank  		   SMALLINT,
	PRIMARY KEY (bioentry_id,dbxref_id)
);

CREATE INDEX dblink_dbx  ON bioentry_dbxref(dbxref_id);


CREATE TABLE reference (
  	reference_id       INT NOT NULL IDENTITY,
	dbxref_id	   INT,
  	location 	   LONGVARCHAR NOT NULL,
  	title    	   LONGVARCHAR,
  	authors  	   LONGVARCHAR NOT NULL,
  	crc	   	   VARCHAR(32),
        PRIMARY KEY (reference_id),
	UNIQUE (dbxref_id),
	UNIQUE (crc)
);


CREATE TABLE bioentry_reference (
  	bioentry_id 	INT NOT NULL,
  	reference_id 	INT NOT NULL,
  	start_pos	INT,
  	end_pos	  	INT,
  	rank  		SMALLINT DEFAULT 0 NOT NULL,
  	PRIMARY KEY(bioentry_id,reference_id,rank)
);

CREATE INDEX bioentryref_ref ON bioentry_reference(reference_id);


-- We use the table name "anncomment" instead of "comment" (which is a reserved word in some RDBMS)
CREATE TABLE anncomment (
  	comment_id  	INT NOT NULL IDENTITY,
  	bioentry_id    	INT NOT NULL,
  	comment_text   	LONGVARCHAR NOT NULL,
  	rank   		SMALLINT DEFAULT 0 NOT NULL,
        PRIMARY KEY (comment_id),
  	UNIQUE(bioentry_id, rank)
);


CREATE TABLE bioentry_qualifier_value (
	bioentry_id   		INT NOT NULL,
   	term_id  		INT NOT NULL,
   	value         		LONGVARCHAR,
	rank			INT DEFAULT 0 NOT NULL,
	UNIQUE (bioentry_id,term_id,rank)
);

CREATE INDEX bioentryqual_trm ON bioentry_qualifier_value(term_id);


CREATE TABLE seqfeature (
   	seqfeature_id 		INT NOT NULL IDENTITY,
   	bioentry_id   		INT NOT NULL,
   	type_term_id		INT NOT NULL,
   	source_term_id  	INT NOT NULL,
	display_name		VARCHAR(64),
   	rank 			SMALLINT DEFAULT 0 NOT NULL,
        PRIMARY KEY (seqfeature_id),
	UNIQUE (bioentry_id,type_term_id,source_term_id,rank)
);

CREATE INDEX seqfeature_trm  ON seqfeature(type_term_id);
CREATE INDEX seqfeature_fsrc ON seqfeature(source_term_id);


CREATE TABLE seqfeature_relationship (
        seqfeature_relationship_id INT NOT NULL IDENTITY,
   	object_seqfeature_id	INT NOT NULL,
   	subject_seqfeature_id 	INT NOT NULL,
   	term_id 	INT NOT NULL,
   	rank 			INT,
        PRIMARY KEY (seqfeature_relationship_id),
	UNIQUE (object_seqfeature_id,subject_seqfeature_id,term_id)
);

CREATE INDEX seqfeaturerel_trm   ON seqfeature_relationship(term_id);
CREATE INDEX seqfeaturerel_child ON seqfeature_relationship(subject_seqfeature_id);


CREATE TABLE seqfeature_path (
   	object_seqfeature_id	INT NOT NULL,
   	subject_seqfeature_id 	INT NOT NULL,
   	term_id 		INT NOT NULL,
	distance	     	INT,
	UNIQUE (object_seqfeature_id,subject_seqfeature_id,term_id,distance)
);

CREATE INDEX seqfeaturepath_trm   ON seqfeature_path(term_id);
CREATE INDEX seqfeaturepath_child ON seqfeature_path(subject_seqfeature_id);


CREATE TABLE seqfeature_qualifier_value (
	seqfeature_id 		INT NOT NULL,
   	term_id 		INT NOT NULL,
   	rank 			SMALLINT DEFAULT 0 NOT NULL,
   	value  			LONGVARCHAR NOT NULL,
   	PRIMARY KEY (seqfeature_id,term_id,rank)
);

CREATE INDEX seqfeaturequal_trm ON seqfeature_qualifier_value(term_id);
   

CREATE TABLE seqfeature_dbxref ( 
       	seqfeature_id      INT NOT NULL,
       	dbxref_id          INT NOT NULL,
  	rank  		   SMALLINT,
	PRIMARY KEY (seqfeature_id,dbxref_id)
);

CREATE INDEX feadblink_dbx  ON seqfeature_dbxref(dbxref_id);


CREATE TABLE location (
	location_id		INT NOT NULL IDENTITY,
   	seqfeature_id		INT NOT NULL,
	dbxref_id		INT,
	term_id			INT,
   	start_pos              	INT,
   	end_pos                	INT,
   	strand             	TINYINT DEFAULT 0 NOT NULL,
   	rank          		SMALLINT DEFAULT 0 NOT NULL,
        PRIMARY KEY (location_id),
   	UNIQUE (seqfeature_id, rank)
);

CREATE INDEX seqfeatureloc_start ON location(start_pos, end_pos);
CREATE INDEX seqfeatureloc_dbx   ON location(dbxref_id);
CREATE INDEX seqfeatureloc_trm   ON location(term_id);


CREATE TABLE location_qualifier_value (
	location_id		INT NOT NULL,
   	term_id 		INT NOT NULL,
   	value  			VARCHAR(255) NOT NULL,
   	int_value 		INT,
	PRIMARY KEY (location_id,term_id)
);

CREATE INDEX locationqual_trm ON location_qualifier_value(term_id);

--
-- Create the foreign key constraints
--

-- ontology term

ALTER TABLE term ADD CONSTRAINT FKont_term
	FOREIGN KEY (ontology_id) REFERENCES ontology(ontology_id)
	ON DELETE CASCADE;

-- term synonyms

ALTER TABLE term_synonym ADD CONSTRAINT FKterm_syn
	FOREIGN KEY (term_id) REFERENCES term(term_id)
	ON DELETE CASCADE;

-- term_dbxref

ALTER TABLE term_dbxref ADD CONSTRAINT FKdbxref_trmdbxref
       	FOREIGN KEY (dbxref_id) REFERENCES dbxref(dbxref_id)
	ON DELETE CASCADE;
ALTER TABLE term_dbxref ADD CONSTRAINT FKterm_trmdbxref
      FOREIGN KEY (term_id) REFERENCES term(term_id)
	ON DELETE CASCADE;

-- term_relationship

ALTER TABLE term_relationship ADD CONSTRAINT FKtrmsubject_trmrel
	FOREIGN KEY (subject_term_id) REFERENCES term(term_id)
	ON DELETE CASCADE;
ALTER TABLE term_relationship ADD CONSTRAINT FKtrmpredicate_trmrel
       	FOREIGN KEY (predicate_term_id) REFERENCES term(term_id)
	ON DELETE CASCADE;
ALTER TABLE term_relationship ADD CONSTRAINT FKtrmobject_trmrel
       	FOREIGN KEY (object_term_id) REFERENCES term(term_id)
	ON DELETE CASCADE;
ALTER TABLE term_relationship ADD CONSTRAINT FKterm_trmrel
       	FOREIGN KEY (ontology_id) REFERENCES ontology(ontology_id)
	ON DELETE CASCADE;

-- term_relationship_term

ALTER TABLE term_relationship_term ADD CONSTRAINT FKtrmrel_trmreltrm
      FOREIGN KEY (term_relationship_id) REFERENCES term_relationship(term_relationship_id)
      ON DELETE CASCADE ;
ALTER TABLE term_relationship_term ADD CONSTRAINT FKtrm_trmreltrm
      FOREIGN KEY (term_id) REFERENCES term(term_id)
      ON DELETE CASCADE ;

-- term_path

ALTER TABLE term_path ADD CONSTRAINT FKtrmsubject_trmpath
	FOREIGN KEY (subject_term_id) REFERENCES term(term_id)
	ON DELETE CASCADE;
ALTER TABLE term_path ADD CONSTRAINT FKtrmpredicate_trmpath
       	FOREIGN KEY (predicate_term_id) REFERENCES term(term_id)
	ON DELETE CASCADE;
ALTER TABLE term_path ADD CONSTRAINT FKtrmobject_trmpath
       	FOREIGN KEY (object_term_id) REFERENCES term(term_id)
	ON DELETE CASCADE;
ALTER TABLE term_path ADD CONSTRAINT FKontology_trmpath
       	FOREIGN KEY (ontology_id) REFERENCES ontology(ontology_id)
	ON DELETE CASCADE;

-- taxon, taxon_name

-- unfortunately, we can't constrain parent_taxon_id as it is violated
-- occasionally by the downloads available from NCBI
-- ALTER TABLE taxon ADD CONSTRAINT FKtaxon_taxon
--         FOREIGN KEY (parent_taxon_id) REFERENCES taxon(taxon_id);
ALTER TABLE taxon_name ADD CONSTRAINT FKtaxon_taxonname
        FOREIGN KEY (taxon_id) REFERENCES taxon(taxon_id)
        ON DELETE CASCADE;

-- bioentry

ALTER TABLE bioentry ADD CONSTRAINT FKtaxon_bioentry
	FOREIGN KEY (taxon_id) REFERENCES taxon(taxon_id)
        ON DELETE CASCADE;
ALTER TABLE bioentry ADD CONSTRAINT FKbiodatabase_bioentry
	FOREIGN KEY (biodatabase_id) REFERENCES biodatabase(biodatabase_id)
        ON DELETE CASCADE;

-- bioentry_relationship

ALTER TABLE bioentry_relationship ADD CONSTRAINT FKterm_bioentryrel
	FOREIGN KEY (term_id) REFERENCES term(term_id)
        ON DELETE CASCADE;
ALTER TABLE bioentry_relationship ADD CONSTRAINT FKparentent_bioentryrel
	FOREIGN KEY (object_bioentry_id) REFERENCES bioentry(bioentry_id)
	ON DELETE CASCADE;
ALTER TABLE bioentry_relationship ADD CONSTRAINT FKchildent_bioentryrel
	FOREIGN KEY (subject_bioentry_id) REFERENCES bioentry(bioentry_id)
	ON DELETE CASCADE;

-- bioentry_path

ALTER TABLE bioentry_path ADD CONSTRAINT FKterm_bioentrypath
	FOREIGN KEY (term_id) REFERENCES term(term_id)
        ON DELETE CASCADE;
ALTER TABLE bioentry_path ADD CONSTRAINT FKparentent_bioentrypath
	FOREIGN KEY (object_bioentry_id) REFERENCES bioentry(bioentry_id)
	ON DELETE CASCADE;
ALTER TABLE bioentry_path ADD CONSTRAINT FKchildent_bioentrypath
	FOREIGN KEY (subject_bioentry_id) REFERENCES bioentry(bioentry_id)
	ON DELETE CASCADE;

-- biosequence

ALTER TABLE biosequence ADD CONSTRAINT FKbioentry_bioseq
	FOREIGN KEY (bioentry_id) REFERENCES bioentry(bioentry_id)
	ON DELETE CASCADE;

-- comment

ALTER TABLE anncomment ADD CONSTRAINT FKbioentry_comment
	FOREIGN KEY(bioentry_id) REFERENCES bioentry(bioentry_id)
	ON DELETE CASCADE;

-- bioentry_dbxref

ALTER TABLE bioentry_dbxref ADD CONSTRAINT FKbioentry_dblink
        FOREIGN KEY (bioentry_id) REFERENCES bioentry(bioentry_id)
	ON DELETE CASCADE;
ALTER TABLE bioentry_dbxref ADD CONSTRAINT FKdbxref_dblink
       	FOREIGN KEY (dbxref_id) REFERENCES dbxref(dbxref_id)
	ON DELETE CASCADE;

-- dbxref_qualifier_value

ALTER TABLE dbxref_qualifier_value ADD CONSTRAINT FKtrm_dbxrefqual
	FOREIGN KEY (term_id) REFERENCES term(term_id)
	ON DELETE CASCADE;
ALTER TABLE dbxref_qualifier_value ADD CONSTRAINT FKdbxref_dbxrefqual
	FOREIGN KEY (dbxref_id) REFERENCES dbxref(dbxref_id)
	ON DELETE CASCADE;

-- bioentry_reference

ALTER TABLE bioentry_reference ADD CONSTRAINT FKbioentry_entryref
	FOREIGN KEY (bioentry_id) REFERENCES bioentry(bioentry_id)
	ON DELETE CASCADE;
ALTER TABLE bioentry_reference ADD CONSTRAINT FKreference_entryref
	FOREIGN KEY (reference_id) REFERENCES reference(reference_id)
	ON DELETE CASCADE;

-- bioentry_qualifier_value

ALTER TABLE bioentry_qualifier_value ADD CONSTRAINT FKbioentry_entqual
	FOREIGN KEY (bioentry_id) REFERENCES bioentry(bioentry_id)
	ON DELETE CASCADE;
ALTER TABLE bioentry_qualifier_value ADD CONSTRAINT FKterm_entqual
	FOREIGN KEY (term_id) REFERENCES term(term_id)
	ON DELETE CASCADE;

-- reference 
ALTER TABLE reference ADD CONSTRAINT FKdbxref_reference
      FOREIGN KEY ( dbxref_id ) REFERENCES dbxref ( dbxref_id ) ;

-- seqfeature

ALTER TABLE seqfeature ADD CONSTRAINT FKterm_seqfeature
	FOREIGN KEY (type_term_id) REFERENCES term(term_id)
	ON DELETE CASCADE;
ALTER TABLE seqfeature ADD CONSTRAINT FKsourceterm_seqfeature
	FOREIGN KEY (source_term_id) REFERENCES term(term_id)
	ON DELETE CASCADE;
ALTER TABLE seqfeature ADD CONSTRAINT FKbioentry_seqfeature
	FOREIGN KEY (bioentry_id) REFERENCES bioentry(bioentry_id)
	ON DELETE CASCADE;

-- seqfeature_relationship

ALTER TABLE seqfeature_relationship ADD CONSTRAINT FKterm_seqfeatrel
	FOREIGN KEY (term_id) REFERENCES term(term_id)
	ON DELETE CASCADE;
ALTER TABLE seqfeature_relationship ADD CONSTRAINT FKparentfeat_seqfeatrel
	FOREIGN KEY (object_seqfeature_id) REFERENCES seqfeature(seqfeature_id)
	ON DELETE CASCADE;
ALTER TABLE seqfeature_relationship ADD CONSTRAINT FKchildfeat_seqfeatrel
	FOREIGN KEY (subject_seqfeature_id) REFERENCES seqfeature(seqfeature_id)
	ON DELETE CASCADE;

-- seqfeature_path

ALTER TABLE seqfeature_path ADD CONSTRAINT FKterm_seqfeatpath
	FOREIGN KEY (term_id) REFERENCES term(term_id)
	ON DELETE CASCADE;
ALTER TABLE seqfeature_path ADD CONSTRAINT FKparentfeat_seqfeatpath
	FOREIGN KEY (object_seqfeature_id) REFERENCES seqfeature(seqfeature_id)
	ON DELETE CASCADE;
ALTER TABLE seqfeature_path ADD CONSTRAINT FKchildfeat_seqfeatpath
	FOREIGN KEY (subject_seqfeature_id) REFERENCES seqfeature(seqfeature_id)
	ON DELETE CASCADE;

-- seqfeature_qualifier_value
ALTER TABLE seqfeature_qualifier_value ADD CONSTRAINT FKterm_featqual
	FOREIGN KEY (term_id) REFERENCES term(term_id)
	ON DELETE CASCADE;
ALTER TABLE seqfeature_qualifier_value ADD CONSTRAINT FKseqfeature_featqual
	FOREIGN KEY (seqfeature_id) REFERENCES seqfeature(seqfeature_id)
	ON DELETE CASCADE;

-- seqfeature_dbxref

ALTER TABLE seqfeature_dbxref ADD CONSTRAINT FKseqfeature_feadblink
        FOREIGN KEY (seqfeature_id) REFERENCES seqfeature(seqfeature_id)
	ON DELETE CASCADE;
ALTER TABLE seqfeature_dbxref ADD CONSTRAINT FKdbxref_feadblink
       	FOREIGN KEY (dbxref_id) REFERENCES dbxref(dbxref_id)
	ON DELETE CASCADE;

-- location

ALTER TABLE location ADD CONSTRAINT FKseqfeature_location
	FOREIGN KEY (seqfeature_id) REFERENCES seqfeature(seqfeature_id)
	ON DELETE CASCADE;
ALTER TABLE location ADD CONSTRAINT FKdbxref_location
	FOREIGN KEY (dbxref_id) REFERENCES dbxref(dbxref_id)
	ON DELETE CASCADE;
ALTER TABLE location ADD CONSTRAINT FKterm_featloc
	FOREIGN KEY (term_id) REFERENCES term(term_id)
	ON DELETE CASCADE;

-- location_qualifier_value

ALTER TABLE location_qualifier_value ADD CONSTRAINT FKfeatloc_locqual
	FOREIGN KEY (location_id) REFERENCES location(location_id)
	ON DELETE CASCADE;
ALTER TABLE location_qualifier_value ADD CONSTRAINT FKterm_locqual
	FOREIGN KEY (term_id) REFERENCES term(term_id)
	ON DELETE CASCADE;

