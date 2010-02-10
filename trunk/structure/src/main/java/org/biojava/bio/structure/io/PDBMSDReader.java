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
 * Created on 20.09.2004
 * @author Andreas Prlic
 *
 */


package org.biojava.bio.structure.io;

import java.io.IOException;
import java.sql.Connection;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;

import javax.sql.DataSource;

import org.biojava.bio.structure.AminoAcidImpl;
import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.AtomImpl;
import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.ChainImpl;
import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.HetatomImpl;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureImpl;
import org.biojava.utils.JDBCPooledDataSource;


/** retreive a structure from an installation of the MSD - search
    database ( http://www.ebi.ac.uk/msd/index.html ,
    http://www.ebi.ac.uk/msd-srv/docs/dbdoc/ ) */
public class PDBMSDReader 
    implements StructureIO {
  
    protected DataSource dataSource ;
   
    
    public PDBMSDReader(){
	super() ;
	dataSource = null ;
    }
    
    /** open a database conenction to the MSD  search database
	@param dbDriver specified which JDBC driver to use e.g. 
	@param dbUrl connection string 
	@param dbUsername the username in the db
	@param dbPassword password to use
    @throws SQLException
    */
    public void setDBConnection( String dbDriver,String dbUrl,String dbUsername, String dbPassword)
	throws SQLException
    {
	try {
	    dataSource = JDBCPooledDataSource.getDataSource(dbDriver, dbUrl, dbUsername, dbPassword);
	}catch (Exception e){
	    e.printStackTrace();
	    throw new SQLException("could not get DataSource from JDBCPooledDataSource");
	}
	
    }
    
    /** Get a structure by providing a PDB code.
     * expects connections parameters to be set a system properties.
     * 
     * 
     * @param pdbId  a String specifying the id value (PDB code)
     * @return a Structure object, or null if no structure with matching PDB code has been found
     * @throws IOException ...
     */

    public Structure getStructureById(String pdbId) throws IOException {

	// open connection to database
	Connection conn ;
	Structure structure = new StructureImpl();
	try {
	    conn = dataSource.getConnection();

	    // create Structure object from MSD database query.
	} catch (SQLException e) {
	    e.printStackTrace();
	    throw new IOException ("could not open database connection") ;
	}

	try {
	    PreparedStatement ps = conn.prepareStatement("select residue_serial, serial, residue_id,code_3_letter, CHEM_ATOM_NAME_PDB_LS, chem_atom_name,        element_symbol,alt_code,  CHAIN_CODE_1_LETTER, chain_code, chain_pdb_code,  OCCUPANCY, PDB_CHARGE,PDB_Group, x,y,z, CHEM_COMP_CODE, CHEM_ATOM_ID,U_ISO_OR_EQUIV,  RESIDUE_PDB_SEQ, RESIDUE_PDB_INSERT_CODE from Atom_Data where accession_code = ? order by serial") ;
	    ps.setString(1,pdbId);
	    System.out.println(ps);
	    ResultSet row = ps.executeQuery();
	    int prev_serial = -9999 ;
	    String prevChain = "" ;
	    String prevType  = "" ; // was it amino acid or hetatom before ?
	    Group g =null;
	    Chain current_chain = new ChainImpl();
	    while (row.next()) {
	    
		//int residue_serial         = row.getInt(1);
		int serial                 = row.getInt(2);
		//String residue_id          = row.getString(3);
		String code_3_letter       = row.getString(4);
		String CHEM_ATOM_NAME_LS   = row.getString(5);
		String chem_atom_name      = row.getString(6);
		//String element_symbol      = row.getString(7);
		//String alt_code            = row.getString(8);
		//String chain_Code_1_Letter = row.getString(9);
		//String chain_code          = row.getString(10) ;
		String chain_pdb_code      = row.getString(11) ;
		double OCCUPANCY           = row.getDouble(12);
		//String PDB_CHARGE          = row.getString(13);
		String pdb_Group           = row.getString(14);
		double x                   = row.getDouble(15);
		double y                   = row.getDouble(16);
		double z                   = row.getDouble(17);
		//String ligand_code 	   = row.getString(18);
		//String chem_atom_id  	   = row.getString(19);
		//double U_ISO_OR_EQUIV      = row.getDouble(20);
		int residue_pdb_seq        = row.getInt(21);
		String insertionCode       = row.getString(22);
		//System.out.println(U_ISO_OR_EQUIV);
	    
	    
		//String str = "" ;
	    
	    


		if ( 
		    ( prev_serial != residue_pdb_seq ) || 
		    ( ! prevType.equals(pdb_Group)      )
		    ) {
		    if ( prev_serial  != -9999) {
			current_chain.addGroup(g);
		    }
		    if ( pdb_Group.equals("A") ) {
			g = new AminoAcidImpl();
		    } else {
			g = new HetatomImpl();
		    }
		}

		if ( insertionCode == null ) 
		    insertionCode = "";

	    	if ( chain_pdb_code == null ) 
		    chain_pdb_code = " " ;
		
		if (! prevChain.equals(chain_pdb_code)) 
		    if ( prevChain != "" ) {
			structure.addChain(current_chain) ;
			current_chain = new ChainImpl();
		    }
		
		
		
		if ( CHEM_ATOM_NAME_LS == null) 
		    CHEM_ATOM_NAME_LS = "";
		
		
		//select residue_serial, chain_code_1_letter,RESIDUE_PDB_CODE, RESIDUE_PDB_INSERT_CODE, SP_PRIMARY_ID, DSC_TYPE, SP_SERIAL from swiss_prot_mapping where accession_code = '5pti' order by ASSEMBLY_SERIAL ;
		//str += residue_serial+" " + serial + " " +code_3_letter+" " +CHEM_ATOM_NAME_LS  + chem_atom_name +" " + chem_atom_id + " " +chain_Code_1_Letter + " " + OCCUPANCY +" " +PDB_CHARGE+" " + x+" "+ y+" "+z + " " +ligand_code ;
		//System.out.println(str);
	    

		Atom a = new AtomImpl();
		a.setX(x);
		a.setY(y);
		a.setZ(z);
		a.setName(chem_atom_name);
		a.setPDBserial(serial);
		//a.setAltLoc(alt_code);
		//a.setTempFactor(U_ISO_OR_EQUIV);
		a.setOccupancy(OCCUPANCY);
		String fname = CHEM_ATOM_NAME_LS  + chem_atom_name ;
		//String fullname = fname + (4-fname.length()) * " ";
		String fullname = fname ;
		for ( int i = 0 ; i < 4-fname.length(); i++ ) {
		    fullname +=" " ;
		}
		//System.out.println(">"+fullname+"<");
		a.setFullName(fullname);
		g.addAtom(a);
		g.setPDBFlag(true);
		g.setPDBName(code_3_letter);
		g.setPDBCode(residue_pdb_seq+insertionCode);
		//System.out.println(a);
		//System.out.println("chainname "+chain_pdb_code+ " " +chain_Code_1_Letter+ " " + chain_code);
		current_chain.setName(chain_pdb_code);
		prev_serial = residue_pdb_seq;
		prevChain   = chain_pdb_code  ;
		prevType    = pdb_Group ;
	    }
	    //System.out.println("adding!" + g + current_chain);
	    if ( g != null) {
		current_chain.addGroup(g);
		structure.addChain(current_chain) ;
	    } else {
		// no structure found 
		return null ;
	    }
	} catch ( Exception e ){	    
	    e.printStackTrace();
	  
	    try {
		conn.close();
	    } catch ( SQLException es) {
		es.printStackTrace();
	    }
	    return null;
	}
    
	try {
	    conn.close();
	} catch ( SQLException es) {
	    es.printStackTrace();
	}

	return structure;
    
       
    }

}


