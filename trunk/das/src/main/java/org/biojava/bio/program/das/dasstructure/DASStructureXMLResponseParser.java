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
 * Created on 19.03.2004
 * @author Andreas Prlic
 *
 */
package org.biojava.bio.program.das.dasstructure ;

import java.util.ArrayList;
import java.util.HashMap;

import org.biojava.bio.seq.ProteinTools;
import org.biojava.bio.seq.io.SymbolTokenization;
import org.biojava.bio.structure.AminoAcidImpl;
import org.biojava.bio.structure.AtomImpl;
import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.ChainImpl;
import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.HetatomImpl;
import org.biojava.bio.structure.NucleotideImpl;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureImpl;
import org.biojava.bio.structure.io.PDBFileReader;
import org.biojava.bio.structure.io.PDBParseException;
import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.Symbol;
import org.xml.sax.Attributes;
import org.xml.sax.helpers.DefaultHandler;


/** a class to Parse the XML response of a DAS structure service. 
 * returns a Structure object
 */
public class DASStructureXMLResponseParser  extends DefaultHandler{

    /**
     * 
     */
    StructureImpl structure;

    ArrayList     current_model  ;
    ChainImpl     current_chain  ;
    Group         current_group  ;

    // for conversion 3code 1code
    SymbolTokenization threeLetter ;
    SymbolTokenization oneLetter ;

    public DASStructureXMLResponseParser() {
	super();
	structure = new StructureImpl() ;
	current_model = new ArrayList();
	current_chain = null           ;
	current_group = null           ;


	// for conversion
	Alphabet alpha_prot = ProteinTools.getAlphabet();

	try {
	    threeLetter = alpha_prot.getTokenization("name");
	    oneLetter  = alpha_prot.getTokenization("token");
	} catch (Exception e) {
	    e.printStackTrace() ;
	}

    }
    
    /**
     * returns the Structure object.
     * @return a Structure object
     */
    public Structure get_structure() {
	return structure ;
    }

 
	

    private void OBJECThandler(Attributes atts){
	// mandatory attribute
	//System.out.println("objectHandler");
	String id = atts.getValue("dbAccessionId") ;
	String objectVersion = atts.getValue("objectVersion") ;

	HashMap header = new HashMap() ;
	header.put("idCode",id);
	header.put("modDate",objectVersion);	
	structure.setHeader(header);
	//structure.setName(id);
	structure.setPDBCode(id);
    }
    


    private void CHAINhandler(Attributes atts){

	String model_nr = atts.getValue("model");
	if ( model_nr != null ) {
	    int modnr = Integer.parseInt(model_nr);
	    //System.out.println("modnr: " + modnr) ;
	    // if DAS server sets modelnr -> a NMR structure ...
	    structure.setNmr(true);

	    if ( modnr -1 > structure.nrModels() ) {
		// a new model starts !
		structure.addModel(current_model);
		current_model = new ArrayList();	    
	    }
	}
	
	String chain_id = atts.getValue("id");
	String spCode   = atts.getValue("SwissprotId");
	//System.out.println("new chain "+chain_id);
	if ( current_chain ==null) {
	    // first chain
	    current_chain = new ChainImpl();
	    current_chain.setName(chain_id);
	    current_chain.setSwissprotId(spCode);
	    return ;
	}
	


	// not the first chain
	// paranoic: check if we had it already ..
	
	Chain testchain = isKnownChain(chain_id);
	
	if (testchain != null) {
	    //System.out.println("already known..."+ chain_id + testchain.getLength());
	    current_chain = (ChainImpl)testchain ;
	    
	} else {
	    current_chain = new ChainImpl();
	    current_chain.setName(chain_id);
	    current_chain.setSwissprotId(spCode);
	    
	}
	
    }

    /**
       parse a line like this: 
       <ATOM x="xCoord" y="yCoord" z="zCoord" atomName="atomname" atomID="atomID" occupancy="occupancy" tempFactor="tempFactor"/>
    */
    
    private void ATOMhandler(Attributes atts){
	String atomID   = atts.getValue("atomID")  ;
	String atomName = atts.getValue("atomName");
	String X        = atts.getValue("x");
	String Y        = atts.getValue("y");
	String Z        = atts.getValue("z");

	// optional occupancy, tempFactor 
	// todo: implement ...
	
	AtomImpl atom = new AtomImpl() ;
	int pdbnumber = Integer.parseInt(atomID);
	atom.setPDBserial(pdbnumber);
	atom.setFullName(atomName);
	atom.setName(atomName.trim());
	double x = Double.parseDouble(X);
	double y = Double.parseDouble(Y);
	double z = Double.parseDouble(Z);
	double[] coords = new double[3];
	
	coords[0] = x ;
	coords[1] = y ;
	coords[2] = z ;
	atom.setCoords(coords);
	current_group.addAtom(atom);

    }
    /* initiale new group, either Hetatom or AminoAcid */
    private Group getNewGroup(String type,String name) {
	Group group;

	if ( type.equals("amino")){
	    AminoAcidImpl aa = new AminoAcidImpl() ;

	    Character aminoCode1 = null;
	    try{
		aminoCode1 = convert_3code_1code(name);
	    } 
	    catch (IllegalSymbolException e) {
		aminoCode1 = new Character('x');
	    }
	    aa.setAminoType(aminoCode1);

	    group = aa ;
	} else if ( type.equals("nucleotide")) {
	    // it is a nucleotidee
	    NucleotideImpl nu = new NucleotideImpl();
	    group = nu;
	} else {
	    group = new HetatomImpl();
	}
	
	try{
	    group.setPDBName(name);
	} catch (PDBParseException e) {
	    
	}
	return group;
    }



    private void GROUPhandler(Attributes atts){

	String name    = atts.getValue("name"   );
	String type    = atts.getValue("type"   );
	String groupid = atts.getValue("groupID");
	
	current_group = getNewGroup(type,name);		   
	current_group.setPDBCode(groupid);
    }
    
	
    public void startElement (String uri, String name, String qName, Attributes atts){
	//System.out.println("new element uri: >"+uri+"< name:>"+name+"< qname:>" +qName+"<");
	if      (qName.equals("object")) OBJECThandler(atts) ;
	else if (qName.equals("chain") ) CHAINhandler (atts) ;
	else if (qName.equals("atom")  ) ATOMhandler  (atts) ;
	else if (qName.equals("group") ) GROUPhandler (atts) ;
					
	    
	    
	
				
    }
	
    public void startDocument() {
	//System.out.println("start document");
	
    }
	
    public void endDocument ()	{
	//features.add(feature);		
	structure.addModel(current_model);
	
    }

    public void endElement(String uri, String name, String qName) {
	//System.out.println("end >"+name+"<");

	if ( qName.equals("group")){
	    current_chain.addGroup(current_group);
	}

	if ( qName.equals("chain")) {
	    // check if chain is already known ...

	    Chain ch = isKnownChain(current_chain.getName());
	    if ( ch == null) {
		current_model.add(current_chain);
	    }
	}

	
    }
    
    public void characters (char ch[], int start, int length){
	//System.out.println("characters");
	//for (int i = start; i < start + length; i++) {
	    
	    //characterdata += ch[i];
	//}
	
    }
    
    /** test if the chain is already known (is in current_model
     * ArrayList) and if yes, returns the chain 
     * @see PDBFileReader
     */
    private Chain isKnownChain(String chainID){
	Chain testchain = null;
	Chain retchain  = null;
	//System.out.println("isKnownCHain: >"+chainID+"< current_chains:"+current_model.size());

	for (int i = 0; i< current_model.size();i++){
	    testchain = (Chain) current_model.get(i);
	    //System.out.println("comparing chainID >"+chainID+"< against testchain " + i+" >" +testchain.getName()+"<");
	    if (chainID.equals(testchain.getName())) {
		//System.out.println("chain "+ chainID+" already known ...");
		retchain = testchain;
		break ;
	    }
	}
	//if (retchain == null) {
	//  System.out.println("unknownCHain!");
	//}
	return retchain;
    }
  

    /** convert three character amino acid codes into single character
     *  e.g. convert CYS to C 
     * @see PDBFileReader
     */
    
    private Character convert_3code_1code(String code3) 
	throws IllegalSymbolException
    {
	Symbol sym   =  threeLetter.parseToken(code3) ;
	String code1 =  oneLetter.tokenizeSymbol(sym);
	
	return new Character(code1.charAt(0)) ;
		
    }


}
