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
 */
package org.biojava.nbio.structure.io.mmcif;

import org.biojava.nbio.structure.io.FileParsingParameters;
import org.biojava.nbio.structure.io.mmcif.chem.MetalBondDistance;
import org.biojava.nbio.structure.io.mmcif.model.*;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Created by andreas on 6/9/16.
 */
public class MetalBondConsumer implements MMcifConsumer{


    Map<String,List<MetalBondDistance>> definitions = new HashMap<>();

    @Override
    public void documentStart() {
        definitions.clear();
    }

    @Override
    public void documentEnd() {

        // minimize memory consumption

        for  (List<MetalBondDistance> d : definitions.values()){
            ArrayList<MetalBondDistance> a = (ArrayList<MetalBondDistance>)d;

            a.trimToSize();
        }

    }

    @Override
    public void newAtomSite(AtomSite atom) {

    }

    @Override
    public void newEntity(Entity entity) {

    }

    @Override
    public void newEntityPoly(EntityPoly entityPoly) {

    }

    @Override
    public void newEntityPolySeq(EntityPolySeq epolseq) {

    }

    @Override
    public void newStructAsym(StructAsym sasym) {

    }

    @Override
    public void setStruct(Struct struct) {

    }

    @Override
    public void newDatabasePDBrev(DatabasePDBrev dbrev) {

    }

    @Override
    public void newDatabasePDBrevRecord(DatabasePdbrevRecord dbrev) {

    }

    @Override
    public void newDatabasePDBremark(DatabasePDBremark remark) {

    }

    @Override
    public void newExptl(Exptl exptl) {

    }

    @Override
    public void newCell(Cell cell) {

    }

    @Override
    public void newSymmetry(Symmetry symmetry) {

    }

    @Override
    public void newStructNcsOper(StructNcsOper sNcsOper) {

    }
    
    @Override 
    public void newAtomSites(AtomSites atomSites) {
    	
    }

    @Override
    public void newStructRef(StructRef sref) {

    }

    @Override
    public void newStructRefSeq(StructRefSeq sref) {

    }

    @Override
    public void newStructRefSeqDif(StructRefSeqDif sref) {

    }

    @Override
    public void newStructSite(StructSite sref) {

    }

    @Override
    public void newStructSiteGen(StructSiteGen sref) {

    }

    @Override
    public void newPdbxPolySeqScheme(PdbxPolySeqScheme ppss) {

    }

    @Override
    public void newPdbxNonPolyScheme(PdbxNonPolyScheme ppss) {

    }

    @Override
    public void newPdbxEntityNonPoly(PdbxEntityNonPoly pen) {

    }

    @Override
    public void newStructKeywords(StructKeywords kw) {

    }

    @Override
    public void newRefine(Refine r) {

    }

    @Override
    public void newChemComp(ChemComp c) {

    }

    @Override
    public void newChemCompDescriptor(ChemCompDescriptor ccd) {

    }

    @Override
    public void newPdbxStructOperList(PdbxStructOperList structOper) {

    }

    @Override
    public void newPdbxStrucAssembly(PdbxStructAssembly strucAssembly) {

    }

    @Override
    public void newPdbxStrucAssemblyGen(PdbxStructAssemblyGen strucAssembly) {

    }

    @Override
    public void newChemCompAtom(ChemCompAtom atom) {

    }

    @Override
    public void newPdbxChemCompIndentifier(PdbxChemCompIdentifier id) {

    }

    @Override
    public void newChemCompBond(ChemCompBond bond) {

    }

    @Override
    public void newPdbxChemCompDescriptor(PdbxChemCompDescriptor desc) {

    }

    @Override
    public void newEntitySrcGen(EntitySrcGen entitySrcGen) {

    }

    @Override
    public void newEntitySrcNat(EntitySrcNat entitySrcNat) {

    }

    @Override
    public void newEntitySrcSyn(EntitySrcSyn entitySrcSyn) {

    }

    @Override
    public void newStructConn(StructConn structConn) {

    }

    @Override
    public void newAuditAuthor(AuditAuthor aa) {

    }

    @Override
    public void newGenericData(String category, List<String> loopFields, List<String> lineData) {

        MetalBondDistance d = new MetalBondDistance();

        d.setAtomType1(lineData.get(0));
        d.setAtomType2(lineData.get(1));
        d.setLowerLimit(Float.parseFloat(lineData.get(2)));
        d.setUpperLimit(Float.parseFloat(lineData.get(3)));

        List<MetalBondDistance> defs = definitions.get(d.getAtomType1());

        if ( defs == null){
            defs = new ArrayList<>();
            definitions.put(d.getAtomType1(),defs);
        }

        defs.add(d);

    }

    @Override
    public void setFileParsingParameters(FileParsingParameters params) {

    }

    @Override
    public FileParsingParameters getFileParsingParameters() {
        return null;
    }

    public Map<String,List<MetalBondDistance>> getDefinitions(){
        return definitions;
    }

	@Override
	public void newPdbxAuditRevisionHistory(PdbxAuditRevisionHistory history) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void newPdbxDatabaseStatus(PdbxDatabaseStatus status) {
		// TODO Auto-generated method stub
		
	}
}
