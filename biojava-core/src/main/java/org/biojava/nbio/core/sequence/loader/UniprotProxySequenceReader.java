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
 * Created on 01-21-2010
 *
 * @auther Scooter Willis
 *
 */
package org.biojava.nbio.core.sequence.loader;

import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.AccessionID;
import org.biojava.nbio.core.sequence.DataSource;
import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.Strand;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompoundSet;
import org.biojava.nbio.core.sequence.features.DBReferenceInfo;
import org.biojava.nbio.core.sequence.features.DatabaseReferenceInterface;
import org.biojava.nbio.core.sequence.features.FeaturesKeyWordInterface;
import org.biojava.nbio.core.sequence.storage.SequenceAsStringHelper;
import org.biojava.nbio.core.sequence.template.*;
import org.biojava.nbio.core.util.Equals;
import org.biojava.nbio.core.util.XMLHelper;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.xml.sax.SAXException;

import javax.xml.parsers.ParserConfigurationException;
import javax.xml.xpath.XPathExpressionException;
import java.io.*;
import java.net.HttpURLConnection;
import java.net.URL;
import java.rmi.RemoteException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.regex.Pattern;

/**
 *
 * Pass in a Uniprot ID and this ProxySequenceReader when passed to a ProteinSequence will get the sequence data and other data elements
 * associated with the ProteinSequence by Uniprot. This is an example of how to map external databases of proteins and features to the BioJava3
 * ProteinSequence.
 * Important to call @see setUniprotDirectoryCache to allow caching of XML files so they don't need to be reloaded each time. Does
 * not manage cache.
 * @param <C>
 */
public class UniprotProxySequenceReader<C extends Compound> implements ProxySequenceReader<C>, FeaturesKeyWordInterface, DatabaseReferenceInterface {

	private final static Logger logger = LoggerFactory.getLogger(UniprotProxySequenceReader.class);

	/*
	 * Taken from http://www.uniprot.org/help/accession_numbers
	 */
	private static final String SPID_PATTERN = "[OPQ][0-9][A-Z0-9]{3}[0-9]";
	private static final String TREMBLID_PATTERN = "[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}";
	public static final Pattern UP_AC_PATTERN = Pattern.compile("(" + SPID_PATTERN + "|" + TREMBLID_PATTERN + ")");

	public static final String DEFAULT_UNIPROT_BASE_URL = "https://www.uniprot.org";

	private static String uniprotbaseURL = DEFAULT_UNIPROT_BASE_URL;
	private static String uniprotDirectoryCache = null;
	private String sequence;
	private CompoundSet<C> compoundSet;
	private List<C> parsedCompounds = new ArrayList<C>();
	Document uniprotDoc;

	/**
	 * The UniProt id is used to retrieve the UniProt XML which is then parsed as a DOM object
	 * so we know everything about the protein. If an error occurs throw an exception. We could
	 * have a bad uniprot id or network error
	 * @param accession
	 * @param compoundSet
	 * @throws CompoundNotFoundException
	 * @throws IOException if problems while reading the UniProt XML
	 */
	public UniprotProxySequenceReader(String accession, CompoundSet<C> compoundSet) throws CompoundNotFoundException, IOException {
		if (!UP_AC_PATTERN.matcher(accession.toUpperCase()).matches()) {
			throw new IllegalArgumentException("Accession provided " + accession + " doesn't comply with the uniprot acession pattern.");
		}
		setCompoundSet(compoundSet);
		uniprotDoc = this.getUniprotXML(accession);
		String seq = this.getSequence(uniprotDoc);
		setContents(seq);
	}

	/**
	 * The xml is passed in as a DOM object so we know everything about the protein.
	 *  If an error occurs throw an exception. We could have a bad uniprot id
	 * @param document
	 * @param compoundSet
	 * @throws CompoundNotFoundException
	 */
	public UniprotProxySequenceReader(Document document, CompoundSet<C> compoundSet) throws CompoundNotFoundException {
		setCompoundSet(compoundSet);
		uniprotDoc = document;
		String seq = this.getSequence(uniprotDoc);
		setContents(seq);
	}
	/**
	 * The passed in xml is parsed as a DOM object so we know everything about the protein.
	 *  If an error occurs throw an exception. We could have a bad uniprot id
	 * @param xml
	 * @param compoundSet
	 * @return UniprotProxySequenceReader
	 * @throws Exception
	 */
	public static <C extends Compound> UniprotProxySequenceReader<C> parseUniprotXMLString(String xml, CompoundSet<C> compoundSet) {
		try {
			Document document = XMLHelper.inputStreamToDocument(new ByteArrayInputStream(xml.getBytes()));
			return new UniprotProxySequenceReader<C>(document, compoundSet);
		} catch (Exception e) {
			logger.error("Exception on xml parse of: {}", xml);
		}
		return null;
	}

	@Override
	public void setCompoundSet(CompoundSet<C> compoundSet) {
		this.compoundSet = compoundSet;
	}

	/**
	 * Once the sequence is retrieved set the contents and make sure everything this is valid
	 * @param sequence
	 * @throws CompoundNotFoundException
	 */
	@Override
	public void setContents(String sequence) throws CompoundNotFoundException {
		// Horrendously inefficient - pretty much the way the old BJ did things.
		// TODO Should be optimised.
		this.sequence = sequence;
		this.parsedCompounds.clear();
		for (int i = 0; i < sequence.length();) {
			String compoundStr = null;
			C compound = null;
			for (int compoundStrLength = 1; compound == null && compoundStrLength <= compoundSet.getMaxSingleCompoundStringLength(); compoundStrLength++) {
				compoundStr = sequence.substring(i, i + compoundStrLength);
				compound = compoundSet.getCompoundForString(compoundStr);
			}
			if (compound == null) {
				throw new CompoundNotFoundException("Compound "+compoundStr+" not found");
			} else {
				i += compoundStr.length();
			}
			this.parsedCompounds.add(compound);
		}
	}

	/**
	 * The sequence length
	 * @return
	 */
	@Override
	public int getLength() {
		return this.parsedCompounds.size();
	}

	/**
	 *
	 * @param position
	 * @return
	 */
	@Override
	public C getCompoundAt(int position) {
		return this.parsedCompounds.get(position - 1);
	}

	/**
	 *
	 * @param compound
	 * @return
	 */
	@Override
	public int getIndexOf(C compound) {
		return this.parsedCompounds.indexOf(compound) + 1;
	}

	/**
	 *
	 * @param compound
	 * @return
	 */
	@Override
	public int getLastIndexOf(C compound) {
		return this.parsedCompounds.lastIndexOf(compound) + 1;
	}

	/**
	 *
	 * @return
	 */
	@Override
	public String toString() {
		return getSequenceAsString();
	}

	/**
	 *
	 * @return
	 */
	@Override
	public String getSequenceAsString() {
		return sequence;
	}

	/**
	 *
	 * @return
	 */
	@Override
	public List<C> getAsList() {
		return this.parsedCompounds;
	}

	@Override
	public boolean equals(Object o){

		if(! Equals.classEqual(this, o)) {
			return false;
		}

		Sequence<C> other = (Sequence<C>)o;
		if ( other.getCompoundSet() != getCompoundSet())
			return false;

		List<C> rawCompounds = getAsList();
		List<C> otherCompounds = other.getAsList();

		if ( rawCompounds.size() != otherCompounds.size())
			return false;

		for (int i = 0 ; i < rawCompounds.size() ; i++){
			Compound myCompound = rawCompounds.get(i);
			Compound otherCompound = otherCompounds.get(i);
			if ( ! myCompound.equalsIgnoreCase(otherCompound))
				return false;
		}
		return true;
	}

	@Override
	public int hashCode(){
		String s = getSequenceAsString();
		return s.hashCode();
	}

	/**
	 *
	 * @return
	 */
	@Override
	public SequenceView<C> getInverse() {
		return SequenceMixin.inverse(this);
	}

	/**
	 *
	 * @param bioBegin
	 * @param bioEnd
	 * @param strand
	 * @return
	 */
	public String getSequenceAsString(Integer bioBegin, Integer bioEnd, Strand strand) {
		SequenceAsStringHelper<C> sequenceAsStringHelper = new SequenceAsStringHelper<C>();
		return sequenceAsStringHelper.getSequenceAsString(this.parsedCompounds, compoundSet, bioBegin, bioEnd, strand);
	}

	/**
	 *
	 * @param bioBegin
	 * @param bioEnd
	 * @return
	 */
	@Override
	public SequenceView<C> getSubSequence(final Integer bioBegin, final Integer bioEnd) {
		return new SequenceProxyView<C>(UniprotProxySequenceReader.this, bioBegin, bioEnd);
	}

	/**
	 *
	 * @return
	 */
	@Override
	public Iterator<C> iterator() {
		return this.parsedCompounds.iterator();
	}

	/**
	 *
	 * @return
	 */
	@Override
	public CompoundSet<C> getCompoundSet() {
		return compoundSet;
	}

	/**
	 *
	 * @return
	 */
	@Override
	public AccessionID getAccession() {
		AccessionID accessionID = new AccessionID();
		if (uniprotDoc == null) {
			return accessionID;
		}
		try {
			Element uniprotElement = uniprotDoc.getDocumentElement();
			Element entryElement = XMLHelper.selectSingleElement(uniprotElement, "entry");
			Element nameElement = XMLHelper.selectSingleElement(entryElement, "name");
			accessionID = new AccessionID(nameElement.getTextContent(), DataSource.UNIPROT);
		} catch (XPathExpressionException e) {
			logger.error("Exception: ", e);
		}
		return accessionID;
	}

	/**
	 * Pull uniprot accessions associated with this sequence
	 * @return
	 * @throws XPathExpressionException
	 */
	public ArrayList<AccessionID> getAccessions() throws XPathExpressionException {
		ArrayList<AccessionID> accessionList = new ArrayList<AccessionID>();
		if (uniprotDoc == null) {
			return accessionList;
		}
		Element uniprotElement = uniprotDoc.getDocumentElement();
		Element entryElement = XMLHelper.selectSingleElement(uniprotElement, "entry");
		ArrayList<Element> keyWordElementList = XMLHelper.selectElements(entryElement, "accession");
		for (Element element : keyWordElementList) {
			AccessionID accessionID = new AccessionID(element.getTextContent(), DataSource.UNIPROT);
			accessionList.add(accessionID);
		}

		return accessionList;
	}

	/**
	 * Pull uniprot protein aliases associated with this sequence
	 * Provided for backwards compatibility now that we support both
	 * gene and protein aliases via separate methods.
	 * @return
	 * @throws XPathExpressionException
	 */
	public ArrayList<String> getAliases() throws XPathExpressionException {

		return getProteinAliases();
	}
	/**
	 * Pull uniprot protein aliases associated with this sequence
	 * @return
	 * @throws XPathExpressionException
	 */
	public ArrayList<String> getProteinAliases() throws XPathExpressionException {
		ArrayList<String> aliasList = new ArrayList<String>();
		if (uniprotDoc == null) {
			return aliasList;
		}
		Element uniprotElement = uniprotDoc.getDocumentElement();
		Element entryElement = XMLHelper.selectSingleElement(uniprotElement, "entry");
		Element proteinElement = XMLHelper.selectSingleElement(entryElement, "protein");
		ArrayList<Element> keyWordElementList = XMLHelper.selectElements(proteinElement, "alternativeName");
		for (Element element : keyWordElementList) {
			Element fullNameElement = XMLHelper.selectSingleElement(element, "fullName");
			aliasList.add(fullNameElement.getTextContent());
			Element shortNameElement = XMLHelper.selectSingleElement(element, "shortName");
			if(null != shortNameElement) {
				String shortName = shortNameElement.getTextContent();
				if(null != shortName && !shortName.trim().isEmpty()) {
					aliasList.add(shortName);
				}
			}
		}
		keyWordElementList = XMLHelper.selectElements(proteinElement, "recommendedName");
		for (Element element : keyWordElementList) {
			Element fullNameElement = XMLHelper.selectSingleElement(element, "fullName");
			aliasList.add(fullNameElement.getTextContent());
			Element shortNameElement = XMLHelper.selectSingleElement(element, "shortName");
			if(null != shortNameElement) {
				String shortName = shortNameElement.getTextContent();
				if(null != shortName && !shortName.trim().isEmpty()) {
					aliasList.add(shortName);
				}
			}
		}
		Element cdAntigen = XMLHelper.selectSingleElement(proteinElement, "cdAntigenName");
		if(null != cdAntigen) {
			String cdAntigenName = cdAntigen.getTextContent();
			if(null != cdAntigenName && !cdAntigenName.trim().isEmpty()) {
				aliasList.add(cdAntigenName);
			}
		}

		return aliasList;
	}

	/**
	 * Pull uniprot gene aliases associated with this sequence
	 * @return
	 * @throws XPathExpressionException
	 */
	public ArrayList<String> getGeneAliases() throws XPathExpressionException {
		ArrayList<String> aliasList = new ArrayList<String>();
		if (uniprotDoc == null) {
			return aliasList;
		}
		Element uniprotElement = uniprotDoc.getDocumentElement();
		Element entryElement = XMLHelper.selectSingleElement(uniprotElement, "entry");
		ArrayList<Element> proteinElements = XMLHelper.selectElements(entryElement, "gene");
		for(Element proteinElement : proteinElements) {
			ArrayList<Element> keyWordElementList = XMLHelper.selectElements(proteinElement, "name");
			for (Element element : keyWordElementList) {
				aliasList.add(element.getTextContent());
			}
		}
		return aliasList;
	}

	/**
	 *
	 * @param compounds
	 * @return
	 */
	@Override
	public int countCompounds(C... compounds) {
		throw new UnsupportedOperationException("Not supported yet.");
	}

	/**
	 *
	 * @param accession
	 * @return
	 * @throws IOException
	 */
	private Document getUniprotXML(String accession) throws IOException, CompoundNotFoundException {
		StringBuilder sb = new StringBuilder();
		// try in cache
		if (uniprotDirectoryCache != null && uniprotDirectoryCache.length() > 0) {
			sb = fetchFromCache(accession);
		}

		// http://www.uniprot.org/uniprot/?query=SORBIDRAFT_03g027040&format=xml
		if (sb.length() == 0) {
			String uniprotURL = getUniprotbaseURL() + "/uniprot/" + accession.toUpperCase() + ".xml";
			logger.info("Loading: {}", uniprotURL);
			sb = fetchUniprotXML(uniprotURL);

			int index = sb.indexOf("xmlns="); //strip out name space stuff to make it easier on xpath
			if (index != -1) {
				int lastIndex = sb.indexOf(">", index);
				sb.replace(index, lastIndex, "");
			}
			if (uniprotDirectoryCache != null && uniprotDirectoryCache.length() > 0)
				writeCache(sb,accession);
		}

		logger.info("Load complete");
		try {
			//       logger.debug(sb.toString());
			Document document = XMLHelper.inputStreamToDocument(new ByteArrayInputStream(sb.toString().getBytes()));
			return document;
		} catch (SAXException e) {
			logger.error("Exception on xml parse of: {}", sb.toString());
		} catch (ParserConfigurationException e) {
			logger.error("Exception on xml parse of: {}", sb.toString());
		}
		return null;
	}

	private void writeCache(StringBuilder sb, String accession) throws IOException {
		File f = new File(uniprotDirectoryCache + File.separatorChar + accession + ".xml");
		FileWriter fw = new FileWriter(f);
		fw.write(sb.toString());
		fw.close();
	}

	private StringBuilder fetchUniprotXML(String uniprotURL)
			throws IOException, CompoundNotFoundException {

		StringBuilder sb = new StringBuilder();
		URL uniprot = new URL(uniprotURL);
		int attempt = 5;
		List<String> errorCodes = new ArrayList<String>();
		while(attempt > 0) {
			HttpURLConnection uniprotConnection = (HttpURLConnection) uniprot.openConnection();
			uniprotConnection.setRequestProperty("User-Agent", "BioJava");
			uniprotConnection.connect();
			int statusCode = uniprotConnection.getResponseCode();
			if (statusCode == 200) {
				BufferedReader in = new BufferedReader(
						new InputStreamReader(
						uniprotConnection.getInputStream()));
				String inputLine;

				while ((inputLine = in.readLine()) != null) {
					sb.append(inputLine);
				}
				in.close();
				return sb;
			}
			attempt--;
			errorCodes.add(String.valueOf(statusCode));
		}
		throw new RemoteException("Couldn't fetch accession from the url " + uniprotURL + " error codes on 5 attempts are " + errorCodes.toString());
	}

	/**
	 * @param key
	 * @return A string containing the contents of entry specified by key and if not found returns an empty string
	 * @throws FileNotFoundException
	 * @throws IOException
	 */
	private StringBuilder fetchFromCache(String key)
			throws FileNotFoundException, IOException {
		int index;
		File f = new File(uniprotDirectoryCache + File.separatorChar + key + ".xml");
		StringBuilder sb = new StringBuilder();
		if (f.exists()) {
			FileReader fr = new FileReader(f);
			int size = (int) f.length();
			char[] data = new char[size];
			fr.read(data);
			fr.close();
			sb.append(data);
			index = sb.indexOf("xmlns="); //strip out name space stuff to make it easier on xpath
			if (index != -1) {
				int lastIndex = sb.indexOf(">", index);
				sb.replace(index, lastIndex, "");
			}
		}
		return sb;
	}

	/**
	 *
	 * @param uniprotDoc
	 * @return
	 */
	private String getSequence(Document uniprotDoc)  {

		try {
			Element uniprotElement = uniprotDoc.getDocumentElement();
			Element entryElement = XMLHelper.selectSingleElement(uniprotElement, "entry");
			Element sequenceElement = XMLHelper.selectSingleElement(entryElement, "sequence");

			String seqdata = sequenceElement.getTextContent();

			return seqdata;
		} catch (XPathExpressionException e) {
			logger.error("Problems while parsing sequence in UniProt XML: {}. Sequence will be blank.", e.getMessage());
			return "";
		}
	}

	/**
	 * The current UniProt URL to deal with caching issues. www.uniprot.org is load balanced
	 * but you can access pir.uniprot.org directly.
	 * @return the uniprotbaseURL
	 */
	public static String getUniprotbaseURL() {
		return uniprotbaseURL;
	}

	/**
	 * @param aUniprotbaseURL the uniprotbaseURL to set
	 */
	public static void setUniprotbaseURL(String aUniprotbaseURL) {
		uniprotbaseURL = aUniprotbaseURL;
	}

	/**
	 * Local directory cache of XML that can be downloaded
	 * @return the uniprotDirectoryCache
	 */
	public static String getUniprotDirectoryCache() {
		return uniprotDirectoryCache;
	}

	/**
	 * @param aUniprotDirectoryCache the uniprotDirectoryCache to set
	 */
	public static void setUniprotDirectoryCache(String aUniprotDirectoryCache) {
		File f = new File(aUniprotDirectoryCache);
		if (!f.exists()) {
			f.mkdirs();
		}
		uniprotDirectoryCache = aUniprotDirectoryCache;
	}

	public static void main(String[] args) {

		try {
			UniprotProxySequenceReader<AminoAcidCompound> uniprotSequence = new UniprotProxySequenceReader<AminoAcidCompound>("YA745_GIBZE", AminoAcidCompoundSet.getAminoAcidCompoundSet());
			ProteinSequence proteinSequence = new ProteinSequence(uniprotSequence);
			logger.info("Accession: {}", proteinSequence.getAccession().getID());
			logger.info("Sequence: {}", proteinSequence.getSequenceAsString());
		} catch (Exception e) {
			logger.error("Exception: ", e);
		}

	}

	/**
	 * Get the gene name associated with this sequence.
	 * @return
	 */
	public String getGeneName() {
		if (uniprotDoc == null) {
			return "";
		}
		try {
			Element uniprotElement = uniprotDoc.getDocumentElement();
			Element entryElement = XMLHelper.selectSingleElement(uniprotElement, "entry");
			Element geneElement = XMLHelper.selectSingleElement(entryElement, "gene");
			if (geneElement == null) {
				return "";
			}
			Element nameElement = XMLHelper.selectSingleElement(geneElement, "name");
			if (nameElement == null) {
				return "";
			}
			return nameElement.getTextContent();
		} catch (XPathExpressionException e) {
			logger.error("Problems while parsing gene name in UniProt XML: {}. Gene name will be blank.",e.getMessage());
			return "";
		}
	}

	/**
	 * Get the organism name assigned to this sequence
	 * @return
	 */
	public String getOrganismName() {
		if (uniprotDoc == null) {
			return "";
		}
		try {
			Element uniprotElement = uniprotDoc.getDocumentElement();
			Element entryElement = XMLHelper.selectSingleElement(uniprotElement, "entry");
			Element organismElement = XMLHelper.selectSingleElement(entryElement, "organism");
			if (organismElement == null) {
				return "";
			}
			Element nameElement = XMLHelper.selectSingleElement(organismElement, "name");
			if (nameElement == null) {
				return "";
			}
			return nameElement.getTextContent();
		} catch (XPathExpressionException e) {
			logger.error("Problems while parsing organism name in UniProt XML: {}. Organism name will be blank.",e.getMessage());
			return "";
		}

	}

	/**
	 * Pull UniProt key words which is a mixed bag of words associated with this sequence
	 * @return
	 */
	@Override
	public ArrayList<String> getKeyWords() {
		ArrayList<String> keyWordsList = new ArrayList<String>();
		if (uniprotDoc == null) {
			return keyWordsList;
		}
		try {
			Element uniprotElement = uniprotDoc.getDocumentElement();

			Element entryElement = XMLHelper.selectSingleElement(uniprotElement, "entry");
			ArrayList<Element> keyWordElementList = XMLHelper.selectElements(entryElement, "keyword");
			for (Element element : keyWordElementList) {
				keyWordsList.add(element.getTextContent());
			}
		} catch (XPathExpressionException e) {
			logger.error("Problems while parsing keywords in UniProt XML: {}. No keywords will be available.",e.getMessage());
			return new ArrayList<String>();
		}

		return keyWordsList;
	}

	/**
	 * The Uniprot mappings to other database identifiers for this sequence
	 * @return
	 */
	@Override
	public LinkedHashMap<String, ArrayList<DBReferenceInfo>> getDatabaseReferences()  {
		LinkedHashMap<String, ArrayList<DBReferenceInfo>> databaseReferencesHashMap = new LinkedHashMap<String, ArrayList<DBReferenceInfo>>();
		if (uniprotDoc == null) {
			return databaseReferencesHashMap;
		}

		try {
			Element uniprotElement = uniprotDoc.getDocumentElement();
			Element entryElement = XMLHelper.selectSingleElement(uniprotElement, "entry");
			ArrayList<Element> dbreferenceElementList = XMLHelper.selectElements(entryElement, "dbReference");
			for (Element element : dbreferenceElementList) {
				String type = element.getAttribute("type");
				String id = element.getAttribute("id");
				ArrayList<DBReferenceInfo> idlist = databaseReferencesHashMap.get(type);
				if (idlist == null) {
					idlist = new ArrayList<DBReferenceInfo>();
					databaseReferencesHashMap.put(type, idlist);
				}
				DBReferenceInfo dbreferenceInfo = new DBReferenceInfo(type, id);
				ArrayList<Element> propertyElementList = XMLHelper.selectElements(element, "property");
				for (Element propertyElement : propertyElementList) {
					String propertyType = propertyElement.getAttribute("type");
					String propertyValue = propertyElement.getAttribute("value");
					dbreferenceInfo.addProperty(propertyType, propertyValue);
				}

				idlist.add(dbreferenceInfo);
			}
		} catch (XPathExpressionException e) {
			logger.error("Problems while parsing db references in UniProt XML: {}. No db references will be available.",e.getMessage());
			return new LinkedHashMap<String, ArrayList<DBReferenceInfo>>();
		}

		return databaseReferencesHashMap;
	}
}
