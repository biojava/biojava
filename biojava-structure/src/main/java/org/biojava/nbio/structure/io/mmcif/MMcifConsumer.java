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
 * created at Mar 4, 2008
 */
package org.biojava.nbio.structure.io.mmcif;

import org.biojava.nbio.structure.io.FileParsingParameters;
import org.biojava.nbio.structure.io.mmcif.model.*;

import java.util.List;

/** An interface for the events triggered by a MMcifParser.
 * The Consumer listens to the events and builds up the protein structure.
 *
 * @author Andreas Prlic
 *  @since 1.7
 *
 */
public interface MMcifConsumer {
	/** called at start of document
	 *
	 */
    void documentStart();

	/** called at end of document
	 *
	 */
    void documentEnd();


	/** A new AtomSite record has been read. Contains the Atom data
	 *
	 * @param atom
	 */
    void newAtomSite(AtomSite atom);
	void newEntity(Entity entity);
	void newEntityPoly(EntityPoly entityPoly);
	void newEntityPolySeq(EntityPolySeq epolseq);
	void newStructAsym(StructAsym sasym);
	void setStruct(Struct struct);
	void newDatabasePDBrev(DatabasePDBrev dbrev);
	void newDatabasePDBrevRecord(DatabasePdbrevRecord dbrev);
	void newDatabasePDBremark(DatabasePDBremark remark);
	void newExptl(Exptl exptl);
	void newCell(Cell cell);
	void newSymmetry(Symmetry symmetry);
	void newStructNcsOper(StructNcsOper sNcsOper);
	void newAtomSites(AtomSites atomSites);
	void newStructRef(StructRef sref);
	void newStructRefSeq(StructRefSeq sref);
	void newStructRefSeqDif(StructRefSeqDif sref);
	void newStructSite(StructSite sref);
	void newStructSiteGen(StructSiteGen sref);
	void newPdbxAuditRevisionHistory(PdbxAuditRevisionHistory history);
	void newPdbxDatabaseStatus(PdbxDatabaseStatus status);
	void newPdbxPolySeqScheme(PdbxPolySeqScheme ppss);
	void newPdbxNonPolyScheme(PdbxNonPolyScheme ppss);
	void newPdbxEntityNonPoly(PdbxEntityNonPoly pen);
	void newStructKeywords(StructKeywords kw);
	void newRefine(Refine r);
	void newChemComp(ChemComp c);
	void newChemCompDescriptor(ChemCompDescriptor ccd);
	void newPdbxStructOperList(PdbxStructOperList structOper);
	void newPdbxStrucAssembly(PdbxStructAssembly strucAssembly);
	void newPdbxStrucAssemblyGen(PdbxStructAssemblyGen strucAssembly);
	void newChemCompAtom(ChemCompAtom atom);
	void newPdbxChemCompIndentifier(PdbxChemCompIdentifier id);
	void newChemCompBond(ChemCompBond bond);
	void newPdbxChemCompDescriptor(PdbxChemCompDescriptor desc);
	void newEntitySrcGen(EntitySrcGen entitySrcGen);
	void newEntitySrcNat(EntitySrcNat entitySrcNat);
	void newEntitySrcSyn(EntitySrcSyn entitySrcSyn);
	void newStructConn(StructConn structConn);

	/** AuditAuthor contains the info from the PDB-AUTHOR records.
	 *
	 * @param aa
	 */
    void newAuditAuthor(AuditAuthor aa);

	/** This method is called if no particular handler for the provided cif category
	 * has been implemented so far.
	 * @param category The category that is being processed.
	 * @param loopFields the fields of this category.
	 * @param lineData the data that is being provided.
	 */
    void newGenericData(String category, List<String> loopFields, List<String> lineData);

	void setFileParsingParameters(FileParsingParameters params);
	FileParsingParameters getFileParsingParameters();








}
