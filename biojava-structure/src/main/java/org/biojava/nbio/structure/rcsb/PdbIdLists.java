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
package org.biojava.nbio.structure.rcsb;

import java.io.*;
import java.net.URL;
import java.net.URLConnection;
import java.net.URLEncoder;
import java.util.*;

/**
 *  Utility classes for retrieving lists of PDB IDs.
 *
 *  @author  Andreas Prlic
 *  @since 4.2.0
 */
public class PdbIdLists {

	/** get the list of current PDB IDs
	 *
	 * @return list of current PDB IDs
	 * @throws IOException
	 */
	public static Set<String> getCurrentPDBIds() throws IOException {
	   String xml ="<orgPdbQuery>\n" +
			   "    <version>head</version>\n" +
			   "    <queryType>org.pdb.query.simple.HoldingsQuery</queryType>\n" +
			   "    <description>Holdings : All Structures</description>\n" +
			   "    <experimentalMethod>ignore</experimentalMethod>\n" +
			   "    <moleculeType>ignore</moleculeType>\n" +
			   "  </orgPdbQuery>";

		return postQuery(xml);
	}


	/** Get the PDB IDs of all virus structures in the current PDB
	 *
	 * @return list of all virus structures
	 * @throws IOException
	 */
	public static Set<String> getAllViruses() throws IOException{
		String xml = "<orgPdbQuery>\n" +
				"        <version>head</version>\n" +
				"        <queryType>org.pdb.query.simple.EntriesOfEntitiesQuery</queryType>\n" +
				"        <description>Entries of :Oligomeric state Search : Min Number of oligomeric state=PAU\n" +
				"        and\n" +
				"        TaxonomyTree Search for Viruses\n" +
				"                </description>\n" +
				"        <parent><![CDATA[<orgPdbCompositeQuery version=\"1.0\">\n" +
				"        <queryRefinement>\n" +
				"        <queryRefinementLevel>0</queryRefinementLevel>\n" +
				"        <orgPdbQuery>\n" +
				"        <version>head</version>\n" +
				"        <queryType>org.pdb.query.simple.BiolUnitQuery</queryType>\n" +
				"        <description>Oligomeric state Search : Min Number of oligomeric state=PAU </description>\n" +
				"        <oligomeric_statemin>PAU</oligomeric_statemin>\n" +
				"        </orgPdbQuery>\n" +
				"        </queryRefinement>\n" +
				"        <queryRefinement>\n" +
				"        <queryRefinementLevel>1</queryRefinementLevel>\n" +
				"        <conjunctionType>and</conjunctionType>\n" +
				"        <orgPdbQuery>\n" +
				"        <version>head</version>\n" +
				"        <queryType>org.pdb.query.simple.TreeEntityQuery</queryType>\n" +
				"        <description>TaxonomyTree Search for Viruses</description>\n" +
				"        <t>1</t>\n" +
				"        <n>10239</n>\n" +
				"        <nodeDesc>Viruses</nodeDesc>\n" +
				"        </orgPdbQuery>\n" +
				"        </queryRefinement>\n" +
				"        </orgPdbCompositeQuery>]]></parent>\n" +
				"        </orgPdbQuery>";

		return postQuery(xml);
	}


	/** get list of all current NMR structures
	 *
	 * @return list of NMR structures
	 * @throws IOException
	 */
	public static Set<String> getNMRStructures() throws IOException{
		String xml = "<orgPdbCompositeQuery version=\"1.0\">\n" +
				" <queryRefinement>\n" +
				"  <queryRefinementLevel>0</queryRefinementLevel>\n" +
				"  <orgPdbQuery>\n" +
				"    <version>head</version>\n" +
				"    <queryType>org.pdb.query.simple.HoldingsQuery</queryType>\n" +
				"    <description>Holdings : All Structures</description>\n" +
				"    <experimentalMethod>ignore</experimentalMethod>\n" +
				"    <moleculeType>ignore</moleculeType>\n" +
				"  </orgPdbQuery>\n" +
				" </queryRefinement>\n" +
				" <queryRefinement>\n" +
				"  <queryRefinementLevel>1</queryRefinementLevel>\n" +
				"  <conjunctionType>and</conjunctionType>\n" +
				"  <orgPdbQuery>\n" +
				"    <version>head</version>\n" +
				"    <queryType>org.pdb.query.simple.ExpTypeQuery</queryType>\n" +
				"    <description>Experimental Method is SOLUTION NMR</description>\n" +
				"    <mvStructure.expMethod.value>SOLUTION NMR</mvStructure.expMethod.value>\n" +
				"    <mvStructure.expMethod.exclusive>y</mvStructure.expMethod.exclusive>\n" +
				"  </orgPdbQuery>\n" +
				" </queryRefinement>\n" +
				"</orgPdbCompositeQuery>\n";


		return postQuery(xml);
	}


	/** get all PDB IDs of gag-polyproteins
	 *
	 * @return list of PDB IDs
	 * @throws IOException
	 */
	public static Set<String> getGagPolyproteins() throws IOException {
		String xml = "<orgPdbCompositeQuery version=\"1.0\">\n" +
				" <queryRefinement>\n" +
				"  <queryRefinementLevel>0</queryRefinementLevel>\n" +
				"  <orgPdbQuery>\n" +
				"    <version>head</version>\n" +
				"    <queryType>org.pdb.query.simple.HoldingsQuery</queryType>\n" +
				"    <description>Holdings : All Structures</description>\n" +
				"    <experimentalMethod>ignore</experimentalMethod>\n" +
				"    <moleculeType>ignore</moleculeType>\n" +
				"  </orgPdbQuery>\n" +
				" </queryRefinement>\n" +
				" <queryRefinement>\n" +
				"  <queryRefinementLevel>1</queryRefinementLevel>\n" +
				"  <conjunctionType>and</conjunctionType>\n" +
				"  <orgPdbQuery>\n" +
				"    <version>head</version>\n" +
				"    <queryType>org.pdb.query.simple.MacroMoleculeQuery</queryType>\n" +
				"    <description>Molecule : Gag-Pol polyprotein [A1Z651, O12158, P03355, P03366, P03367, P03369, P04584, P04585, P04586, P04587, P04588, P05896, P05897, P05959, P05961, P0C6F2, P12497, P12499, P18042, P19505 ... ]</description>\n" +
				"    <macromoleculeName>A1Z651,O12158,P03355,P03366,P03367,P03369,P04584,P04585,P04586,P04587,P04588,P05896,P05897,P05959,P05961,P0C6F2,P12497,P12499,P18042,P19505,P19560,P20875,P24740,P35963,Q699E2,Q70XD7,Q72547,Q7SMT3,Q7SPG9,Q90VT5</macromoleculeName>\n" +
				"  </orgPdbQuery>\n" +
				" </queryRefinement>\n" +
				"</orgPdbCompositeQuery>";

		return postQuery(xml);
	}

	/** get all Transmembrane proteins
	 *
	 * @return list of PDB IDs
	 * @throws IOException
	 */
	public static Set<String> getTransmembraneProteins() throws IOException {
		String xml = "  <orgPdbQuery>\n" +
				"    <version>head</version>\n" +
				"    <queryType>org.pdb.query.simple.TreeQuery</queryType>\n" +
				"    <description>TransmembraneTree Search for root</description>\n" +
				"    <t>19</t>\n" +
				"    <n>0</n>\n" +
				"    <nodeDesc>root</nodeDesc>\n" +
				"  </orgPdbQuery>";

		return postQuery(xml);
	}

	public static Set<String> getNucleotides() throws IOException{
		String xml ="<orgPdbQuery>\n" +
				"    <version>head</version>\n" +
				"    <queryType>org.pdb.query.simple.ChainTypeQuery</queryType>\n" +
				"    <description>Chain Type: there is not any Protein chain</description>\n" +
				"    <containsProtein>N</containsProtein>\n" +
				"    <containsDna>?</containsDna>\n" +
				"    <containsRna>?</containsRna>\n" +
				"    <containsHybrid>?</containsHybrid>\n" +
				"  </orgPdbQuery>";
		return postQuery(xml);
	}

	public static Set<String>getRibosomes() throws IOException{
		String xml = "<orgPdbQuery>\n" +
				"    <version>head</version>\n" +
				"    <queryType>org.pdb.query.simple.StructureKeywordsQuery</queryType>\n" +
				"    <description>StructureKeywordsQuery: struct_keywords.pdbx_keywords.comparator=contains struct_keywords.pdbx_keywords.value=RIBOSOME </description>\n" +
				"    <struct_keywords.pdbx_keywords.comparator>contains</struct_keywords.pdbx_keywords.comparator>\n" +
				"    <struct_keywords.pdbx_keywords.value>RIBOSOME</struct_keywords.pdbx_keywords.value>\n" +
				"  </orgPdbQuery>";

		return postQuery(xml);
	}

	public static final String SERVICELOCATION="http://www.rcsb.org/pdb/rest/search";


	/** post am XML query (PDB XML query format)  to the RESTful RCSB web service
	 *
	 * @param xml
	 * @return a list of PDB ids.
	 */
	public static Set<String> postQuery(String xml)
			throws IOException{

		//System.out.println(xml);


		URL u = new URL(SERVICELOCATION);


		String encodedXML = URLEncoder.encode(xml,"UTF-8");


		InputStream in =  doPOST(u,encodedXML);

		Set<String> pdbIds = new TreeSet<String>();


		try (BufferedReader rd = new BufferedReader(new InputStreamReader(in))) {

			String line;
			while ((line = rd.readLine()) != null) {

				pdbIds.add(line);

			}
			rd.close();
		}


		return pdbIds;



	}

	/** do a POST to a URL and return the response stream for further processing elsewhere.
	 *
	 *
	 * @param url
	 * @return
	 * @throws IOException
	 */
	public static InputStream doPOST(URL url, String data)

			throws IOException
	{

		// Send data

		URLConnection conn = url.openConnection();

		conn.setDoOutput(true);

		try(OutputStreamWriter wr = new OutputStreamWriter(conn.getOutputStream())) {

			wr.write(data);
			wr.flush();
		}


		// Get the response
		return conn.getInputStream();

	};

	public static void main(String[] args){
		try {
			System.out.println("Current PDB status: " + getCurrentPDBIds().size());
			System.out.println("Virus structures: " + getAllViruses().size());
			System.out.println("NMR structures: " + getNMRStructures().size());
			System.out.println("Gag-polyproteins: " + getGagPolyproteins().size());
			System.out.println("Transmembrane proteins: " + getTransmembraneProteins().size());
			System.out.println("Nucleotide: " + getNucleotides().size());
			System.out.println("Ribosomes: " + getRibosomes().size());
		} catch ( Exception e){
			e.printStackTrace();
		}
	}
}
