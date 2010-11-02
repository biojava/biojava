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

public class UniprotProxySequenceReader<C extends Compound> implements ProxySequenceReader<C>, FeaturesKeyWordInterface, DatabaseReferenceInterface {

    private static final Logger logger = Logger.getLogger(UniprotProxySequenceReader.class.getName());
    private static String uniprotbaseURL = "http://pir.uniprot.org";
    private static String uniprotDirectoryCache = null;


    private String sequence;
    private CompoundSet<C> compoundSet;
    private List<C> parsedCompounds = new ArrayList<C>();
    Document uniprotDoc;

    public UniprotProxySequenceReader(String accession, CompoundSet<C> compoundSet) throws Exception {

        setCompoundSet(compoundSet);
        uniprotDoc = this.getUniprotXML(accession);
        String seq = this.getSequence(uniprotDoc);
        setContents(seq);
    }

    public void setCompoundSet(CompoundSet<C> compoundSet) {
        this.compoundSet = compoundSet;
    }

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

    public int getLength() {
        return this.parsedCompounds.size();
    }

    public C getCompoundAt(int position) {
        return this.parsedCompounds.get(position - 1);
    }

    public int getIndexOf(C compound) {
        return this.parsedCompounds.indexOf(compound) + 1;
    }

    public int getLastIndexOf(C compound) {
        return this.parsedCompounds.lastIndexOf(compound) + 1;
    }

    
    public String toString() {
        return getSequenceAsString();
    }

    public String getSequenceAsString() {
        return sequence;
    }

    public List<C> getAsList() {
        return this.parsedCompounds;
    }

    @Override
    public SequenceView<C> getReverse() {
        return SequenceMixin.reverse(this);
    }

    
    public String getSequenceAsString(Integer bioBegin, Integer bioEnd, Strand strand) {
        SequenceAsStringHelper<C> sequenceAsStringHelper = new SequenceAsStringHelper<C>();
        return sequenceAsStringHelper.getSequenceAsString(this.parsedCompounds, compoundSet, bioBegin, bioEnd, strand);
    }

    public SequenceView<C> getSubSequence(final Integer bioBegin, final Integer bioEnd) {
        return new SequenceProxyView<C>(UniprotProxySequenceReader.this, bioBegin, bioEnd);
    }

    public Iterator<C> iterator() {
        return this.parsedCompounds.iterator();
    }

    public CompoundSet<C> getCompoundSet() {
        return compoundSet;
    }

    
    public AccessionID getAccession() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    
    public int countCompounds(C... compounds) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    private Document getUniprotXML(String accession) throws Exception {
        int index = accession.lastIndexOf(".");
        String key = accession;
        if (index != -1) {
            key = accession.substring(0, index);
        }
        StringBuffer sb = new StringBuffer();
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
        if (sb.length() == 0) {
            String uniprotURL = getUniprotbaseURL() + "/uniprot/" + key + ".xml";
            logger.severe("Loading " + uniprotURL);
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

        logger.severe("Load complete");
        try {
            //       System.out.println(sb.toString());
            Document document = XMLHelper.inputStreamToDocument(new ByteArrayInputStream(sb.toString().getBytes()));
            return document;
        } catch (Exception e) {
            System.out.println("Exception on xml parse of:" + sb.toString());
        }
        return null;
    }

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
        if(f.exists() == false){
            f.mkdirs();
        }
        uniprotDirectoryCache = aUniprotDirectoryCache;
    }

    public static void main(String[] args) {

        try {
            UniprotProxySequenceReader<AminoAcidCompound> uniprotSequence = new UniprotProxySequenceReader<AminoAcidCompound>("YA745_GIBZE", AminoAcidCompoundSet.getAminoAcidCompoundSet());
            ProteinSequence proteinSequence = new ProteinSequence(uniprotSequence);

            System.out.println("Sequence=" + proteinSequence.getSequenceAsString());
        } catch (Exception e) {
            e.printStackTrace();
        }

    }

    
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
