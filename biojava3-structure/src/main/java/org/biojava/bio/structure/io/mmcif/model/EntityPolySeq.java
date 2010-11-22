/*
 *                    BioJava development code
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  If you do not have a copy,
 * see:
 *
 *      http://www.gnu.org/copyleft/lesser.html
 *
 * Copyright for this code is held jointly by the individual
 * authors.  These should be listed in @author doc comments.
 *
 * For more information on the BioJava project and its aims,
 * or to join the biojava-l mailing list, visit the home page
 * at:
 *
 *      http://www.biojava.org/
 *
 * created at Jun 1, 2008
 */
package org.biojava.bio.structure.io.mmcif.model;


/** Container for _entity_poly_seq records
 *
<pre>
  Field Name 	 		mmCIF Data Item
  Section   	  		n.a.
  Serial_No   	  		n.a.
  Strand_ID   			PDB strand ID corresponding to   _entity_poly_seq.entity_id   	**
  Strand_Length   	  	derived
  Residue_Names   	  	_entity_poly_seq.mon_id
** Chemically distinct polymer strands are mapped to mmCIF entities. Two instances or the same polymer molecule in the PDB data file are mapped to a single mmCIF entity (eg. a homodimer). For convenience a table of monomer label correspondences is stored in category   PDBX_POLY_SEQ_SCHEME

</pre>
 * @author Andreas Prlic
 * @since 1.7
 */
public class EntityPolySeq extends AbstractBean{
	String entity_id;
	String num;
	String mon_id;
	String hetero;
	public String getEntity_id() {
		return entity_id;
	}
	public void setEntity_id(String entity_id) {
		this.entity_id = entity_id;
	}
	public String getNum() {
		return num;
	}
	public void setNum(String num) {
		this.num = num;
	}
	public String getMon_id() {
		return mon_id;
	}
	public void setMon_id(String mon_id) {
		this.mon_id = mon_id;
	}
	public String getHetero() {
		return hetero;
	}
	public void setHetero(String hetero) {
		this.hetero = hetero;
	}

}
