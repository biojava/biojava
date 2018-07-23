/*
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
 * Created on 26.04.2004
 * @author Andreas Prlic
 *
 */
package org.biojava.nbio.structure.io;

import java.io.IOException;
import java.text.DateFormat;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.text.SimpleDateFormat;
import java.util.List;
import java.util.Locale;

import org.biojava.nbio.core.util.XMLWriter;
import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Bond;
import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.DBRef;
import org.biojava.nbio.structure.Element;
import org.biojava.nbio.structure.Group;
import org.biojava.nbio.structure.GroupType;
import org.biojava.nbio.structure.PDBHeader;
import org.biojava.nbio.structure.Site;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.io.mmcif.MMCIFFileTools;
import org.biojava.nbio.structure.io.mmcif.SimpleMMcifParser;
import org.biojava.nbio.structure.io.mmcif.model.AtomSite;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;


/** 
 * Methods to convert a structure object into different file formats.
 * @author Andreas Prlic
 * @since 1.4
 */
public class FileConvert {

	private static final Logger logger = LoggerFactory.getLogger(FileConvert.class);



	private Structure structure ;

	private boolean printConnections;

	// Locale should be english, e.g. in DE separator is "," -> PDB files have "." !
	public static DecimalFormat d3 = (DecimalFormat)NumberFormat.getInstance(Locale.US);
	static {
		d3.setMaximumIntegerDigits(4);
		d3.setMinimumFractionDigits(3);
		d3.setMaximumFractionDigits(3);
		d3.setGroupingUsed(false);
	}
	public static DecimalFormat d2 = (DecimalFormat)NumberFormat.getInstance(Locale.US);
	static {
		d2.setMaximumIntegerDigits(3);
		d2.setMinimumFractionDigits(2);
		d2.setMaximumFractionDigits(2);
		d2.setGroupingUsed(false);
	}

	private static final String newline = System.getProperty("line.separator");

	/**
	 * Constructs a FileConvert object.
	 *
	 * @param struc  a Structure object
	 */
	public FileConvert(Structure struc) {
		structure = struc ;
		printConnections = true;
	}

	/**
	 * Returns if the Connections should be added
	 * default is true;
	 * @return if the printConnections flag is set
	 */
	public boolean doPrintConnections() {
		return printConnections;
	}

	/** enable/disable printing of connections
	 * connections are sometimes buggy in PDB files
	 * so there are some cases where one might turn this off.
	 * @param printConnections
	 */
	public void setPrintConnections(boolean printConnections) {
		this.printConnections = printConnections;
	}

	/** 
	 * Prints the connections in PDB style
	 *
	 * Rewritten since 5.0 to use {@link Bond}s
	 * Will produce strictly one CONECT record per bond (won't group several bonds in one line)
	 */
	private String printPDBConnections(){

		StringBuilder str = new StringBuilder();

		for (Chain c:structure.getChains()) {
			for (Group g:c.getAtomGroups()) {
				for (Atom a:g.getAtoms()) {
					if (a.getBonds()!=null) {
						for (Bond b:a.getBonds()) {				//7890123456789012345678901234567890123456789012345678901234567890		
							str.append(String.format("CONECT%5d%5d                                                                "+newline, b.getAtomA().getPDBserial(), b.getAtomB().getPDBserial()));
						}
					}
				}
			}
		}
		
		return str.toString();
	}

	/** Convert a structure into a PDB file.
	 * @return a String representing a PDB file.
	 */
	public String toPDB() {


		StringBuffer str = new StringBuffer();
		//int i = 0 ;



		// TODO: print all the PDB header informaton in PDB style
		// some objects (PDBHeader, Compound) are still missing
		//

		PDBHeader header = structure.getPDBHeader();
		header.toPDB(str);


		//REMARK 800
		if (!structure.getSites().isEmpty()) {
			str.append("REMARK 800                                                                      ").append(newline);
			str.append("REMARK 800 SITE                                                                 ").append(newline);
			for (Site site : structure.getSites()) {
				site.remark800toPDB(str);
			}
		}
		//DBREF
		for (DBRef dbref : structure.getDBRefs()){
			dbref.toPDB(str);
			str.append(newline);
		}
		//SSBOND
		List<SSBondImpl> ssbonds = SSBondImpl.getSsBondListFromBondList(structure.getSSBonds());
		for (SSBondImpl ssbond : ssbonds){
			ssbond.toPDB(str);
			str.append(newline);
		}
		//SITE
		for (Site site : structure.getSites()) {
			try {
				site.toPDB(str);
			} catch (Exception e){
				e.printStackTrace();
			}
		}

		//
		// print the atom records
		//

		// do for all models
		int nrModels = structure.nrModels() ;
		if ( structure.isNmr()) {
			str.append("EXPDTA    NMR, "+ nrModels+" STRUCTURES"+newline) ;
		}
		for (int m = 0 ; m < nrModels ; m++) {
			
			
			if ( nrModels>1 ) {
				str.append("MODEL      " + (m+1)+ newline);
			}
			
			List<Chain> polyChains = structure.getPolyChains(m);
			List<Chain> nonPolyChains = structure.getNonPolyChains(m);
			List<Chain> waterChains = structure.getWaterChains(m);
			
			for (Chain chain : polyChains) {

				// do for all groups
				int nrGroups = chain.getAtomLength();
				for ( int h=0; h<nrGroups;h++){

					Group g= chain.getAtomGroup(h);

					toPDB(g,str);

				}
				// End any polymeric chain with a "TER" record
				if (nrGroups > 0) str.append(String.format("%-80s","TER")).append(newline);

			}
			
			boolean nonPolyGroupsExist = false;
			for (Chain chain : nonPolyChains) {

				// do for all groups
				int nrGroups = chain.getAtomLength();
				for ( int h=0; h<nrGroups;h++){

					Group g= chain.getAtomGroup(h);

					toPDB(g,str);
					
					nonPolyGroupsExist = true;
				}

			}
			if (nonPolyGroupsExist) str.append(String.format("%-80s","TER")).append(newline);;

			boolean waterGroupsExist = false;
			for (Chain chain : waterChains) {

				// do for all groups
				int nrGroups = chain.getAtomLength();
				for ( int h=0; h<nrGroups;h++){

					Group g= chain.getAtomGroup(h);

					toPDB(g,str);
					
 					waterGroupsExist = true;
				}

			}
			if (waterGroupsExist) str.append(String.format("%-80s","TER")).append(newline);;


			if ( nrModels>1) {
				str.append(String.format("%-80s","ENDMDL")).append(newline);
			}



		}

		if ( doPrintConnections() )
			str.append(printPDBConnections());

		return str.toString() ;
	}

	private static void toPDB(Group g, StringBuffer str) {
		// iterate over all atoms ...
		// format output ...
		int groupsize  = g.size();

		for ( int atompos = 0 ; atompos < groupsize; atompos++) {
			Atom a = null ;

			a = g.getAtom(atompos);
			if ( a == null)
				continue ;

			toPDB(a, str);


			//line = record + serial + " " + fullname +altLoc
			//+ leftResName + " " + chainID + resseq
			//+ "   " + x+y+z
			//+ occupancy + tempfactor;
			//str.append(line + newline);
			//System.out.println(line);
		}
		if ( g.hasAltLoc()){
			for (Group alt : g.getAltLocs() ) {
				toPDB(alt,str);
			}
		}

	}

	/** Prints the content of an Atom object as a PDB formatted line.
	 *
	 * @param a
	 * @return
	 */
	public static String toPDB(Atom a){
		StringBuffer w = new StringBuffer();

		toPDB(a,w);

		return w.toString();

	}

	public static String toPDB(Atom a, String chainId) {
		StringBuffer w = new StringBuffer();

		toPDB(a,w, chainId);

		return w.toString();
	}


	/**
	 * Convert a Chain object to PDB representation
	 *
	 * @param chain
	 * @return
	 */
	public static String toPDB(Chain chain){
		StringBuffer w = new StringBuffer();
		int nrGroups = chain.getAtomLength();

		for ( int h=0; h<nrGroups;h++){

			Group g= chain.getAtomGroup(h);


			toPDB(g,w);


		}

		return w.toString();
	}

	/**
	 * Convert a Group object to PDB representation
	 *
	 * @param g
	 * @return
	 */
	public static String toPDB(Group g){
		StringBuffer w = new StringBuffer();
		toPDB(g,w);
		return w.toString();
	}

	/**
	 * Print ATOM record in the following syntax
	 * <pre>
	 * ATOM      1  N   ASP A  15     110.964  24.941  59.191  1.00 83.44           N
	 *
	 * COLUMNS        DATA TYPE       FIELD         DEFINITION
	 * ---------------------------------------------------------------------------------
	 * 1 -  6        Record name     "ATOM  "
	 * 7 - 11        Integer         serial        Atom serial number.
	 * 13 - 16        Atom            name          Atom name.
	 * 17             Character       altLoc        Alternate location indicator.
	 * 18 - 20        Residue name    resName       Residue name.
	 * 22             Character       chainID       Chain identifier.
	 * 23 - 26        Integer         resSeq        Residue sequence number.
	 * 27             AChar           iCode         Code for insertion of residues.
	 * 31 - 38        Real(8.3)       x             Orthogonal coordinates for X in
	 * Angstroms.
	 * 39 - 46        Real(8.3)       y             Orthogonal coordinates for Y in
	 * Angstroms.
	 * 47 - 54        Real(8.3)       z             Orthogonal coordinates for Z in
	 * Angstroms.
	 * 55 - 60        Real(6.2)       occupancy     Occupancy.
	 * 61 - 66        Real(6.2)       tempFactor    Temperature factor.
	 * 73 - 76        LString(4)      segID         Segment identifier, left-justified.
	 * 77 - 78        LString(2)      element       Element symbol, right-justified.
	 * 79 - 80        LString(2)      charge        Charge on the atom.
	 * </pre>
	 * @param a
	 * @param str
	 * @param chainID the chain ID that the Atom will have in the output string
	 */
	public static void toPDB(Atom a, StringBuffer str, String chainID) {

		Group g = a.getGroup();

		GroupType type = g.getType() ;

		String record = "" ;
		if ( type.equals(GroupType.HETATM) ) {
			record = "HETATM";
		} else {
			record = "ATOM  ";
		}


		// format output ...
		String resName = g.getPDBName();
		String pdbcode = g.getResidueNumber().toString();


		int    seri       = a.getPDBserial()        ;
		String serial     = String.format("%5d",seri);
		String fullName   = formatAtomName(a);

		Character  altLoc = a.getAltLoc();		
		if ( altLoc == null)
			altLoc = ' ';
		
		String resseq = "" ;
		if ( hasInsertionCode(pdbcode) )
			resseq     = String.format("%5s",pdbcode);
		else
			resseq     = String.format("%4s",pdbcode)+" ";

		String x          = String.format("%8s",d3.format(a.getX()));
		String y          = String.format("%8s",d3.format(a.getY()));
		String z          = String.format("%8s",d3.format(a.getZ()));
		String occupancy  = String.format("%6s",d2.format(a.getOccupancy())) ;
		String tempfactor = String.format("%6s",d2.format(a.getTempFactor()));


		String leftResName = String.format("%3s",resName);

		StringBuffer s = new StringBuffer();
		s.append(record);
		s.append(serial);
		s.append(" ");
		s.append(fullName);
		s.append(altLoc);
		s.append(leftResName);
		s.append(" ");
		s.append(chainID);
		s.append(resseq);
		s.append("   ");
		s.append(x);
		s.append(y);
		s.append(z);
		s.append(occupancy);
		s.append(tempfactor);

		Element e = a.getElement();

		String eString = e.toString().toUpperCase();

		if ( e.equals(Element.R)) {
			eString = "X";
		}
		str.append(String.format("%-76s%2s", s.toString(),eString));
		str.append(newline);

	}

	public static void toPDB(Atom a, StringBuffer str) {
		toPDB(a,str,a.getGroup().getChain().getName());
	}


	/** test if pdbserial has an insertion code */
	private static boolean hasInsertionCode(String pdbserial) {
		try {
			Integer.parseInt(pdbserial) ;
		} catch (NumberFormatException e) {
			return true ;
		}
		return false ;
	}


	/**
	 * Convert a protein Structure to a DAS Structure XML response .
	 * Since 5.0, bond (CONECT records) information is not supported anymore.
	 * @param xw  a XMLWriter object
	 * @throws IOException ...
	 *
	 */
	public void toDASStructure(XMLWriter xw)
			throws IOException
	{

		/*xmlns="http://www.sanger.ac.uk/xml/das/2004/06/17/dasalignment.xsd" xmlns:align="http://www.sanger.ac.uk/xml/das/2004/06/17/alignment.xsd" xmlns:xsd="http://www.w3.org/2001/XMLSchema-instance" xsd:schemaLocation="http://www.sanger.ac.uk/xml/das/2004/06/17/dasalignment.xsd http://www.sanger.ac.uk/xml/das//2004/06/17/dasalignment.xsd"*/

		if ( structure == null){
			System.err.println("can not convert structure null");
			return;
		}

		PDBHeader header = structure.getPDBHeader();

		xw.openTag("object");
		xw.attribute("dbAccessionId",structure.getPDBCode());
		xw.attribute("intObjectId"  ,structure.getPDBCode());
		// missing modification date
		DateFormat dateFormat = new SimpleDateFormat("dd-MMM-yy",Locale.US);
		String modificationDate = dateFormat.format(header.getModDate());
		xw.attribute("objectVersion",modificationDate);
		xw.attribute("type","protein structure");
		xw.attribute("dbSource","PDB");
		xw.attribute("dbVersion","20070116");
		xw.attribute("dbCoordSys","PDBresnum,Protein Structure");

		// do we need object details ???
		xw.closeTag("object");


		// do for all models
		for (int modelnr = 0;modelnr<structure.nrModels();modelnr++){

			// do for all chains:
			for (int chainnr = 0;chainnr<structure.size(modelnr);chainnr++){
				Chain chain = structure.getChainByIndex(modelnr,chainnr);
				xw.openTag("chain");
				xw.attribute("id",chain.getId());
				xw.attribute("SwissprotId",chain.getSwissprotId() );
				if (structure.nrModels()>1){
					xw.attribute("model",Integer.toString(modelnr+1));
				}

				//do for all groups:
				for (int groupnr =0;
						groupnr<chain.getAtomLength()
						;groupnr++){
					Group gr = chain.getAtomGroup(groupnr);
					xw.openTag("group");
					xw.attribute("name",gr.getPDBName());
					xw.attribute("type",gr.getType().toString());
					xw.attribute("groupID",gr.getResidueNumber().toString());


					// do for all atoms:
					//Atom[] atoms  = gr.getAtoms();
					List<Atom> atoms =  gr.getAtoms();
					for (int atomnr=0;atomnr<atoms.size();atomnr++){
						Atom atom = atoms.get(atomnr);
						xw.openTag("atom");
						xw.attribute("atomID",Integer.toString(atom.getPDBserial()));
						xw.attribute("atomName",formatAtomName(atom));
						xw.attribute("x",Double.toString(atom.getX()));
						xw.attribute("y",Double.toString(atom.getY()));
						xw.attribute("z",Double.toString(atom.getZ()));
						xw.closeTag("atom");
					}
					xw.closeTag("group") ;
				}

				xw.closeTag("chain");
			}
		}


		if ( doPrintConnections() ) {
			// not supported anymore since 5.0
		}
	}

	private static String formatAtomName(Atom a) {

		String fullName = null;
		String name = a.getName();
		Element element = a.getElement();

		// RULES FOR ATOM NAME PADDING: 4 columns in total: 13, 14, 15, 16

		// if length 4: nothing to do
		if (name.length()==4)
			fullName = name;

		// if length 3: they stay at 14
		else if (name.length()==3)
			fullName = " "+name;

		// for length 2 it depends:
		//    carbon, oxygens, nitrogens, phosphorous stay at column 14
		//    elements with 2 letters (e.g. NA, FE) will go to column 13
		else if (name.length()==2) {
			if (element == Element.C || element == Element.N || element == Element.O || element == Element.P || element == Element.S)
				fullName = " "+name+" ";
			else
				fullName = name+"  ";
		}

		// for length 1 (e.g. K but also C, O) they stay in column 14
		else if (name.length()==1)
			fullName = " "+name+"  ";

		//if (fullName.length()!=4)
		//	logger.warn("Atom name "+fullName+"to be written in PDB format does not have length 4. Formatting will be incorrect");

		return fullName;
	}


	public String toMMCIF() {

		StringBuilder str = new StringBuilder();

		str.append(SimpleMMcifParser.MMCIF_TOP_HEADER+"BioJava_mmCIF_file"+newline);

		if (structure.getPDBHeader()!=null && structure.getPDBHeader().getCrystallographicInfo()!=null &&
				structure.getPDBHeader().getCrystallographicInfo().getSpaceGroup()!=null &&
				structure.getPDBHeader().getCrystallographicInfo().getCrystalCell()!=null) {

			str.append(MMCIFFileTools.toMMCIF("_cell",
					MMCIFFileTools.convertCrystalCellToCell(structure.getPDBHeader().getCrystallographicInfo().getCrystalCell())));
			str.append(MMCIFFileTools.toMMCIF("_symmetry",
					MMCIFFileTools.convertSpaceGroupToSymmetry(structure.getPDBHeader().getCrystallographicInfo().getSpaceGroup())));

		}


		str.append(getAtomSiteHeader());

		List<AtomSite> list =  MMCIFFileTools.convertStructureToAtomSites(structure);


		str.append(MMCIFFileTools.toMMCIF(list,AtomSite.class));

		return str.toString();
	}

	public static String toMMCIF(Chain chain, String authId, String asymId, boolean writeHeader) {
		StringBuilder str = new StringBuilder();

		if (writeHeader)
			str.append(getAtomSiteHeader());


		List<AtomSite> list = MMCIFFileTools.convertChainToAtomSites(chain, 1, authId, asymId);

		str.append(MMCIFFileTools.toMMCIF(list,AtomSite.class));
		return str.toString();
	}

	public static String toMMCIF(Chain chain, boolean writeHeader) {
		StringBuilder sb = new StringBuilder();
		sb.append(SimpleMMcifParser.MMCIF_TOP_HEADER+"BioJava_mmCIF_file"+newline);
		sb.append(toMMCIF(chain, chain.getName(), chain.getId(),writeHeader));
		return sb.toString();
	}

	public static String getAtomSiteHeader() {
		String header;
		try {
			header = MMCIFFileTools.toLoopMmCifHeaderString("_atom_site", AtomSite.class.getName());

		} catch (ClassNotFoundException e) {
			logger.error("Class not found, will not have a header for this MMCIF category: "+e.getMessage());
			header = "";
		}

		return header;
	}
}
