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
/**
 * This class stores base URL and some cgi arguments for 
 * accessing web-based sequences in NCBI and Swiss-prot, pubmed articles and locuslinks.
 */

package org.biojava.bio.seq.db;


/**
 * @author Lei Lai
 * @author Thomas Down
 * @author Matthew Pocock
 */
public class FetchURL
{
	String baseURL;
	String db;//database name
	String rettype;//return type
	String retmode;//return mode 
	
	/**

	 * Constructs a fetchURL object based on the database name 

	 * and specified return format of sequence.

	 */
	public FetchURL(String databaseName, String format) 
	{	
		if (databaseName.trim().equalsIgnoreCase("genbank")
			||databaseName.trim().equalsIgnoreCase("nucleotide"))
		{
			db = "nucleotide";	
		//	rettype = format;
		//	retmode = format;
		}
		if (databaseName.trim().equalsIgnoreCase("genpept")
			||databaseName.trim().equalsIgnoreCase("protein"))
		{
			db = "protein";
			rettype = format;
			retmode = format;
		}
		if (databaseName.trim().equalsIgnoreCase("swiss-prot"))
		{
			db="swiss-prot";
		}
		if (databaseName.trim().equalsIgnoreCase("pubmed"))
		{
			db="pubmed";
			rettype = "abstract";
			retmode = format;
		}
		if (databaseName.trim().equalsIgnoreCase("locuslink"))
		{
			db="locuslink";
		}
	}
	
	public String getbaseURL()
	{
		if (db.equalsIgnoreCase("Genbank")||db.equalsIgnoreCase("nucleotide")
			||db.equalsIgnoreCase("Genpept")||db.equalsIgnoreCase("protein")
			||db.equalsIgnoreCase("pubmed"))
			baseURL="http://www.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?";
		else if (db.equalsIgnoreCase("Swiss-prot"))
			baseURL="http://us.expasy.org/cgi-bin/get-sprot-raw.pl?";
		else if (db.equalsIgnoreCase("LocusLink"))
			baseURL="http://www.ncbi.nlm.nih.gov/LocusLink/LocRpt.cgi?";
		
		return baseURL;	
	}
	
	//get the database name
	public String getDB()
	{
		return ("db="+db);	
	}
	
	//get the return format and type
	public String getReturnFormat()
	{
		return ("rettype="+rettype+"&retmode="+retmode);	
	}
}
