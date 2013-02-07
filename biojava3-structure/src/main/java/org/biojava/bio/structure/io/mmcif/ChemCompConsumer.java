package org.biojava.bio.structure.io.mmcif;

import java.util.List;

import org.biojava.bio.structure.io.FileParsingParameters;
import org.biojava.bio.structure.io.mmcif.chem.ResidueType;
import org.biojava.bio.structure.io.mmcif.model.AtomSite;
import org.biojava.bio.structure.io.mmcif.model.AuditAuthor;
import org.biojava.bio.structure.io.mmcif.model.ChemComp;
import org.biojava.bio.structure.io.mmcif.model.ChemCompAtom;
import org.biojava.bio.structure.io.mmcif.model.ChemCompBond;
import org.biojava.bio.structure.io.mmcif.model.ChemCompDescriptor;
import org.biojava.bio.structure.io.mmcif.model.DatabasePDBremark;
import org.biojava.bio.structure.io.mmcif.model.DatabasePDBrev;
import org.biojava.bio.structure.io.mmcif.model.Entity;
import org.biojava.bio.structure.io.mmcif.model.EntityPolySeq;
import org.biojava.bio.structure.io.mmcif.model.Exptl;
import org.biojava.bio.structure.io.mmcif.model.PdbxChemCompDescriptor;
import org.biojava.bio.structure.io.mmcif.model.PdbxChemCompIdentifier;
import org.biojava.bio.structure.io.mmcif.model.PdbxEntityNonPoly;
import org.biojava.bio.structure.io.mmcif.model.PdbxNonPolyScheme;
import org.biojava.bio.structure.io.mmcif.model.PdbxPolySeqScheme;
import org.biojava.bio.structure.io.mmcif.model.PdbxStructAssembly;
import org.biojava.bio.structure.io.mmcif.model.PdbxStructAssemblyGen;
import org.biojava.bio.structure.io.mmcif.model.PdbxStructOperList;
import org.biojava.bio.structure.io.mmcif.model.Refine;
import org.biojava.bio.structure.io.mmcif.model.Struct;
import org.biojava.bio.structure.io.mmcif.model.StructAsym;
import org.biojava.bio.structure.io.mmcif.model.StructKeywords;
import org.biojava.bio.structure.io.mmcif.model.StructRef;
import org.biojava.bio.structure.io.mmcif.model.StructRefSeq;

public class ChemCompConsumer implements MMcifConsumer {

	ChemicalComponentDictionary dictionary;

	String latestChemCompId;
	public ChemCompConsumer(){
		dictionary = new ChemicalComponentDictionary();
	}

	public void documentStart() {


	}

	public ChemicalComponentDictionary getDictionary(){
		return dictionary;
	}

	public void newChemComp(ChemComp c) {
		latestChemCompId = c.getId();
		dictionary.addChemComp(c);		
		if ( c.getResidueType() == ResidueType.nonPolymer)
			return;

		if ( c.getResidueType() == ResidueType.saccharide)
			return;

		if ( c.getResidueType() == ResidueType.dSaccharide)
			return;

		//if ( c.isStandard())
		//	System.out.println(c);
	}

	public void documentEnd() {


	}

	public void newAtomSite(AtomSite atom) {
		// TODO Auto-generated method stub

	}

	public void newDatabasePDBremark(DatabasePDBremark remark) {
		// TODO Auto-generated method stub

	}

	public void newDatabasePDBrev(DatabasePDBrev dbrev) {
		// TODO Auto-generated method stub

	}

	public void newEntity(Entity entity) {
		// TODO Auto-generated method stub

	}

	public void newEntityPolySeq(EntityPolySeq epolseq) {
		// TODO Auto-generated method stub

	}

	public void newExptl(Exptl exptl) {
		// TODO Auto-generated method stub

	}

	public void newPdbxEntityNonPoly(PdbxEntityNonPoly pen) {
		// TODO Auto-generated method stub

	}

	public void newPdbxNonPolyScheme(PdbxNonPolyScheme ppss) {
		// TODO Auto-generated method stub

	}

	public void newPdbxPolySeqScheme(PdbxPolySeqScheme ppss) {
		// TODO Auto-generated method stub

	}

	public void newRefine(Refine r) {
		// TODO Auto-generated method stub

	}

	public void newStructAsym(StructAsym sasym) {
		// TODO Auto-generated method stub

	}

	public void newStructKeywords(StructKeywords kw) {
		// TODO Auto-generated method stub

	}

	public void newStructRef(StructRef sref) {
		// TODO Auto-generated method stub

	}

	public void newStructRefSeq(StructRefSeq sref) {
		// TODO Auto-generated method stub

	}

	public void setStruct(Struct struct) {
		// TODO Auto-generated method stub

	}

	public void newGenericData(String category, List<String> loopFields,
			List<String> lineData) {
		//System.out.println("unhandled category: " + category);

	}


	public void newAuditAuthor(AuditAuthor aa)
	{
		// TODO Auto-generated method stub

	}

	public FileParsingParameters getFileParsingParameters()
	{
		// can be ingored in this case...
		return null;
	}

	public void setFileParsingParameters(FileParsingParameters params)
	{
		// TODO Auto-generated method stub

	}

	public void newChemCompDescriptor(ChemCompDescriptor ccd) {
		ChemComp cc = dictionary.getChemComp(latestChemCompId);
		cc.getDescriptors().add(ccd);

	}

	@Override
	public void newPdbxStructOperList(PdbxStructOperList structOper) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void newPdbxStrucAssembly(PdbxStructAssembly strucAssembly) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void newPdbxStrucAssemblyGen(PdbxStructAssemblyGen strucAssembly) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void newChemCompAtom(ChemCompAtom atom) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void newPdbxChemCompIndentifier(PdbxChemCompIdentifier id) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void newChemCompBond(ChemCompBond bond) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void newPdbxChemCompDescriptor(PdbxChemCompDescriptor desc) {
		// TODO Auto-generated method stub
		
	}

}
