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
package org.biojava3.core.sequence.loader;

import java.io.BufferedReader;
import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.InputStreamReader;
import java.net.URL;
import java.net.URLConnection;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.logging.Logger;
import org.biojava3.core.sequence.AccessionID;

import org.biojava3.core.sequence.template.SequenceProxyView;
import org.biojava3.core.sequence.template.Compound;
import org.biojava3.core.exceptions.CompoundNotFoundError;
import org.biojava3.core.sequence.DataSource;
import org.biojava3.core.sequence.ProteinSequence;
import org.biojava3.core.sequence.Strand;
import org.biojava3.core.sequence.compound.AminoAcidCompound;
import org.biojava3.core.sequence.compound.AminoAcidCompoundSet;
import org.biojava3.core.sequence.features.DBReferenceInfo;
import org.biojava3.core.sequence.features.DatabaseReferenceInterface;
import org.biojava3.core.sequence.features.FeaturesKeyWordInterface;

import org.biojava3.core.sequence.storage.SequenceAsStringHelper;
import org.biojava3.core.sequence.template.CompoundSet;
import org.biojava3.core.sequence.template.ProxySequenceReader;
import org.biojava3.core.sequence.template.SequenceMixin;
import org.biojava3.core.sequence.template.SequenceView;
import org.biojava3.core.util.XMLHelper;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

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

    private static final Logger logger = Logger.getLogger(UniprotProxySequenceReader.class.getName());
    private static String uniprotbaseURL = "http://www.uniprot.org"; //"http://pir.uniprot.org";
    private static String uniprotDirectoryCache = null;
    private String sequence;
    private CompoundSet<C> compoundSet;
    private List<C> parsedCompounds = new ArrayList<C>();
    Document uniprotDoc;

    /**
     * The uniprot id is used to retrieve the uniprot XML which is then parsed as a DOM object
     * so we know everything about the protein. If an error occurs throw an exception. We could
     * have a bad uniprot id or network error
     * @param accession
     * @param compoundSet
     * @throws Exception
     */
    public UniprotProxySequenceReader(String accession, CompoundSet<C> compoundSet) throws Exception {

        setCompoundSet(compoundSet);
        uniprotDoc = this.getUniprotXML(accession);
        String seq = this.getSequence(uniprotDoc);
        setContents(seq);
    }

    public void setCompoundSet(CompoundSet<C> compoundSet) {
        this.compoundSet = compoundSet;
    }

    /**
     * Once the sequence is retrieved set the contents and make sure everything this is valid
     * @param sequence
     */
    public void setContents(String sequence) {
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
                throw new CompoundNotFoundError(compoundStr);
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
    public int getLength() {
        return this.parsedCompounds.size();
    }

    /**
     *
     * @param position
     * @return
     */
    public C getCompoundAt(int position) {
        return this.parsedCompounds.get(position - 1);
    }

    /**
     *
     * @param compound
     * @return
     */
    public int getIndexOf(C compound) {
        return this.parsedCompounds.indexOf(compound) + 1;
    }

    /**
     *
     * @param compound
     * @return
     */
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
    public String getSequenceAsString() {
        return sequence;
    }

    /**
     *
     * @return
     */
    public List<C> getAsList() {
        return this.parsedCompounds;
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
    public SequenceView<C> getSubSequence(final Integer bioBegin, final Integer bioEnd) {
        return new SequenceProxyView<C>(UniprotProxySequenceReader.this, bioBegin, bioEnd);
    }

    /**
     *
     * @return
     */
    public Iterator<C> iterator() {
        return this.parsedCompounds.iterator();
    }

    /**
     *
     * @return
     */
    public CompoundSet<C> getCompoundSet() {
        return compoundSet;
    }

    /**
     *
     * @return
     */
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
        } catch (Exception e) {
            e.printStackTrace();
        }
        return accessionID;
    }

    /**
     *
     * @param compounds
     * @return
     */
    public int countCompounds(C... compounds) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    /**
     *
     * @param accession
     * @return
     * @throws Exception
     */
    private Document getUniprotXML(String accession) throws Exception {
        int index = accession.lastIndexOf(".");
        String key = accession;
        if (index != -1) {
            key = accession.substring(0, index);
        }
        StringBuilder sb = new StringBuilder();
        File f = null;
        if (uniprotDirectoryCache != null && uniprotDirectoryCache.length() > 0) {
            f = new File(uniprotDirectoryCache + File.separatorChar + key + ".xml");
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

        }

        // http://www.uniprot.org/uniprot/?query=SORBIDRAFT_03g027040&format=xml
        if (sb.length() == 0) {
            String uniprotURL = getUniprotbaseURL() + "/uniprot/?query=" + key + "&format=xml";
            logger.info("Loading " + uniprotURL);
            URL uniprot = new URL(uniprotURL);
            URLConnection uniprotConnection = uniprot.openConnection();
            BufferedReader in = new BufferedReader(
                    new InputStreamReader(
                    uniprotConnection.getInputStream()));
            String inputLine;

            while ((inputLine = in.readLine()) != null) {
                sb.append(inputLine);
            }
            in.close();
            index = sb.indexOf("xmlns="); //strip out name space stuff to make it easier on xpath
            if (index != -1) {
                int lastIndex = sb.indexOf(">", index);
                sb.replace(index, lastIndex, "");
            }
            if (f != null) {
                FileWriter fw = new FileWriter(f);
                fw.write(sb.toString());
                fw.close();
            }
        }

        logger.info("Load complete");
        try {
            //       System.out.println(sb.toString());
            Document document = XMLHelper.inputStreamToDocument(new ByteArrayInputStream(sb.toString().getBytes()));
            return document;
        } catch (Exception e) {
            System.out.println("Exception on xml parse of:" + sb.toString());
        }
        return null;
    }

    /**
     *
     * @param uniprotDoc
     * @return
     * @throws Exception
     */
    private String getSequence(Document uniprotDoc) throws Exception {
        Element uniprotElement = uniprotDoc.getDocumentElement();
        Element entryElement = XMLHelper.selectSingleElement(uniprotElement, "entry");
        Element sequenceElement = XMLHelper.selectSingleElement(entryElement, "sequence");

        String seqdata = sequenceElement.getTextContent();

        return seqdata;
    }

    /**
     * The current unirpot URL to deal with caching issues. www.uniprot.org is loaded balanced
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
            System.out.println(proteinSequence.getAccession().getID());
            System.out.println("Sequence=" + proteinSequence.getSequenceAsString());
        } catch (Exception e) {
            e.printStackTrace();
        }

    }

    /**
     * Get the gene name associated with this sequence. 
     * @return
     * @throws Exception
     */
    public String getGeneName() throws Exception {
        if (uniprotDoc == null) {
            return "";
        }
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
    }

    /**
     * Get the organism name assigned to this sequence
     * @return
     * @throws Exception
     */
    public String getOrganismName() throws Exception {
        if (uniprotDoc == null) {
            return "";
        }
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
    }

    /**
     * Pull uniprot key words which is a mixed bag of words associated with this sequence
     * @return
     * @throws Exception
     */
    public ArrayList<String> getKeyWords() throws Exception {
        ArrayList<String> keyWordsList = new ArrayList<String>();
        if (uniprotDoc == null) {
            return keyWordsList;
        }
        Element uniprotElement = uniprotDoc.getDocumentElement();
        Element entryElement = XMLHelper.selectSingleElement(uniprotElement, "entry");
        ArrayList<Element> keyWordElementList = XMLHelper.selectElements(entryElement, "keyword");
        for (Element element : keyWordElementList) {
            keyWordsList.add(element.getTextContent());
        }

        return keyWordsList;
    }

    /**
     * The Uniprot mappings to other database identifiers for this sequence
     * @return
     * @throws Exception
     */
    public LinkedHashMap<String, ArrayList<DBReferenceInfo>> getDatabaseReferences() throws Exception {
        LinkedHashMap<String, ArrayList<DBReferenceInfo>> databaseReferencesHashMap = new LinkedHashMap<String, ArrayList<DBReferenceInfo>>();
        if (uniprotDoc == null) {
            return databaseReferencesHashMap;
        }

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


        return databaseReferencesHashMap;
    }
}
