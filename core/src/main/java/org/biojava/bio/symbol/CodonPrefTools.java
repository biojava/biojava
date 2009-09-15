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

package org.biojava.bio.symbol;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.StringTokenizer;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;

import org.biojava.bio.BioException;
import org.biojava.bio.dist.Count;
import org.biojava.bio.dist.Distribution;
import org.biojava.bio.dist.DistributionTools;
import org.biojava.bio.dist.IndexedCount;
import org.biojava.bio.seq.RNATools;
import org.biojava.bio.seq.io.SymbolTokenization;
import org.biojava.utils.ChangeVetoException;
import org.biojava.utils.ClassTools;
import org.biojava.utils.xml.PrettyXMLWriter;
import org.biojava.utils.xml.XMLWriter;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

/**
 * An utility class for codon preferences
 *
 * @author David Huen
 * @author Mark Schreiber
 * @since 1.3
 */
public class CodonPrefTools
{
    /**
     * constants for model organisms
     */
    static String JUNIT = "jUnit use only!!!!";
    /**
     * Drosophila melanogaster codon preferences
     */
    public static String DROSOPHILA_MELANOGASTER_NUCLEAR = "Drosophila melanogaster";
    /**
     * Homo sapiens codon preferences
     */
    public static String MAN_NUCLEAR = "Homo sapiens";
    /**
     * Mus musculus codon preferences
     */
    public static String MOUSE_NUCLEAR = "Mus musculus";
    /**
     * Rattus norvegicus codon preferences
     */
    public static String RAT_NUCLEAR = "Rattus norvegicus";
    /**
     * Takifugu rubripes codon preferences
     */
    public static String FUGU_NUCLEAR = "Takifugu rubripes";
    /**
     * Caenorhabditis elegans codon preferences
     */
    public static String WORM_NUCLEAR = "Caenorhabditis elegans";
    /**
     * Saccharomyces cerevisiae codon preferences
     */
    public static String CEREVISIAE_NUCLEAR = "Saccharomyces cerevisiae";
    /**
     * Schizosaccharomyces pombe codon preferences
     */
    public static String POMBE_NUCLEAR = "Schizosaccharomyces pombe";
    /**
     * Escherichia coli codon preferences
     */
    public static String ECOLI = "Escherichia coli";

    private static Map prefMap;

    final private static AtomicSymbol [] cutg = new AtomicSymbol[64];

    static {
        prefMap = new HashMap();

        loadCodonPreferences();

        try {
            loadCodonOrder();
        }
        catch (IllegalSymbolException ise) {}
    }

    private static class LoadEverythingSelector implements CodonPrefFilter
    {
        public boolean isRequired(String id) { return true; }
        public void put(CodonPref codonPref)
        {
            prefMap.put(codonPref.getName(), codonPref);
        }
    }

    /**
     * get the specified codon preference.
     */
    public static CodonPref getCodonPreference(String id)
    {
        return (CodonPref) prefMap.get(id);
    }

    private static void loadCodonPreferences()
    {
        try {
            // parse the predefined codon preferences
            InputStream prefStream = ClassTools.getClassLoader(CodonPrefTools.class).getResourceAsStream(
                "org/biojava/bio/symbol/CodonPrefTables.xml"
            );

            CodonPrefFilter select = new LoadEverythingSelector();
            readFromXML(prefStream, select);
        }
        catch (Exception e) { e.printStackTrace(); }
    }

    /**
     * returns an RNA dinucleotide alphabet.
     * Used to represent the non-wobble bases in WobbleDistribution
     */
    public static FiniteAlphabet getDinucleotideAlphabet()
    {
        return (FiniteAlphabet)AlphabetManager.generateCrossProductAlphaFromName("(RNA x RNA)");
    }

    /**
     * write out a specified CodonPref object in XML format.
     */
    public static void writeToXML(CodonPref codonPref, PrintWriter writer)
        throws NullPointerException, IOException, IllegalSymbolException, BioException
    {
        XMLWriter xw = new PrettyXMLWriter(writer);

        dumpToXML(codonPref, xw, true);

        writer.flush();
    }

    /**
     * reads a specified CodonPref from an file.
     * @param name name of organism
     */
    public static CodonPref readFromXML(InputStream prefStream, String name)
        throws BioException
    {
        CodonPrefFilter.ByName filter = new CodonPrefFilter.ByName(name);

        readFromXML(prefStream, filter);

        return filter.getCodonPref();
    }

    public static CodonPref[] readFromXML(InputStream prefStream) throws BioException{
      CodonPrefFilter.AcceptAll filter = new CodonPrefFilter.AcceptAll();
      readFromXML(prefStream, filter);

      List l = filter.getCodonPrefs();
      CodonPref[] cp = new CodonPref[l.size()];
      return (CodonPref[])l.toArray(cp);
    }

    /**
     * read an CodonPref XML stream and handle it with a CodonPrefFilter object.
     */
    public static void readFromXML(InputStream prefStream, CodonPrefFilter filter)
        throws BioException
    {
        try {
            DocumentBuilder parser = DocumentBuilderFactory.newInstance().newDocumentBuilder();
            Document doc = parser.parse(prefStream);

            // get tables for each species
            NodeList children = doc.getDocumentElement().getChildNodes();

            for (int i=0; i<children.getLength(); i++) {
                Node cnode = children.item(i);

                if (!(cnode instanceof Element)) continue;

                Element child = (Element) cnode;

                String name = child.getNodeName();

                // the node must be a CodonPref record
                if (!name.equals("CodonPref")) continue;

                // pick up the id and genetic code
                String codonPrefId = child.getAttribute("id");
                String geneticCodeId = child.getAttribute("geneticCodeId");

                // is this entry one we want?
                if (!filter.isRequired(codonPrefId)) continue;

                // now handle each codon frequency entry
                NodeList freqs = child.getChildNodes();

                // create a Count object for the job
                Count freqCounts = new IndexedCount(RNATools.getCodonAlphabet());

                for (int j=0; j < freqs.getLength(); j++) {
                    // load each entry
                    Node freq = freqs.item(j);

                    if (!(freq instanceof Element)) continue;

                    Element freqElement = (Element) freq;

                    // get attributes
                    String codonString = freqElement.getAttribute("codon");
                    String freqString = freqElement.getAttribute("value");

                    // create codon
                    SymbolList codonSL = RNATools.createRNA(codonString);

                    if (codonSL.length() !=3) throw new BioException("'" + codonString + "' is not a valid codon!");

                    AtomicSymbol codon = (AtomicSymbol) RNATools.getCodonAlphabet().getSymbol(codonSL.toList());

                    // recover frequency value too
                    double freqValue = Double.parseDouble(freqString);
                    freqCounts.increaseCount(codon, freqValue);

                }

                // turn the Counts into a Distribution
                Distribution freqDistribution = DistributionTools.countToDistribution(freqCounts);

                // create a CodonPref object
                CodonPref newCodonPref = new SimpleCodonPref(geneticCodeId, freqDistribution, codonPrefId);

                filter.put(newCodonPref);
            }
        }
        catch (Exception e) {
            throw new BioException(e);
        }
    }

    /**
     * reads in a file in Codon Usage Database format and
     * translate it into our XML format
     * These can be obtained from the
     * <a href="http://www.kazusa.or.jp/codon/">Codon Usage Database</a>.
     * <p>
     * Note that the output assumes that the universal genetic code is
     * used as that is not encoded in the CUD files.  Edit the output appropriately
     * to modify the genetic code if necessary.
     */
    public static void translateCUD(InputStream input, OutputStream output)
        throws IOException
    {
        // create a BufferedReader for the job
        BufferedReader rdr = new BufferedReader(new InputStreamReader(input));

        // create a PrintWriter for the job
        PrintWriter pw = new PrintWriter(output);
        CodonPrefFilter.EverythingToXML filter = new CodonPrefFilter.EverythingToXML(pw);

        // now invoke the CUD reader and stream its output to the XML writer
        readFromCUD(rdr, filter);

        filter.close();
    }


    /**
     * converts a String representation of a codon to its Symbol
     */
    private static AtomicSymbol getCodon(String codonString)
        throws IllegalSymbolException
    {
        return (AtomicSymbol) RNATools.getCodonAlphabet().getSymbol(RNATools.createRNA(codonString).toList());
    }

    private static void loadCodonOrder()
        throws IllegalSymbolException
    {
        cutg[0] = getCodon("cga");
        cutg[1] = getCodon("cgc");
        cutg[2] = getCodon("cgg");
        cutg[3] = getCodon("cgu");

        cutg[4] = getCodon("aga");
        cutg[5] = getCodon("agg");

        cutg[6] = getCodon("cua");
        cutg[7] = getCodon("cuc");
        cutg[8] = getCodon("cug");
        cutg[9] = getCodon("cuu");

        cutg[10] = getCodon("uua");
        cutg[11] = getCodon("uug");

        cutg[12] = getCodon("uca");
        cutg[13] = getCodon("ucc");
        cutg[14] = getCodon("ucg");
        cutg[15] = getCodon("ucu");

        cutg[16] = getCodon("agc");
        cutg[17] = getCodon("agu");

        cutg[18] = getCodon("aca");
        cutg[19] = getCodon("acc");
        cutg[20] = getCodon("acg");
        cutg[21] = getCodon("acu");

        cutg[22] = getCodon("cca");
        cutg[23] = getCodon("ccc");
        cutg[24] = getCodon("ccg");
        cutg[25] = getCodon("ccu");

        cutg[26] = getCodon("gca");
        cutg[27] = getCodon("gcc");
        cutg[28] = getCodon("gcg");
        cutg[29] = getCodon("gcu");

        cutg[30] = getCodon("gga");
        cutg[31] = getCodon("ggc");
        cutg[32] = getCodon("ggg");
        cutg[33] = getCodon("ggu");

        cutg[34] = getCodon("gua");
        cutg[35] = getCodon("guc");
        cutg[36] = getCodon("gug");
        cutg[37] = getCodon("guu");

        cutg[38] = getCodon("aaa");
        cutg[39] = getCodon("aag");

        cutg[40] = getCodon("aac");
        cutg[41] = getCodon("aau");

        cutg[42] = getCodon("caa");
        cutg[43] = getCodon("cag");

        cutg[44] = getCodon("cac");
        cutg[45] = getCodon("cau");

        cutg[46] = getCodon("gaa");
        cutg[47] = getCodon("gag");

        cutg[48] = getCodon("gac");
        cutg[49] = getCodon("gau");

        cutg[50] = getCodon("uac");
        cutg[51] = getCodon("uau");

        cutg[52] = getCodon("ugc");
        cutg[53] = getCodon("ugu");

        cutg[54] = getCodon("uuc");
        cutg[55] = getCodon("uuu");

        cutg[56] = getCodon("aua");
        cutg[57] = getCodon("auc");
        cutg[58] = getCodon("auu");

        cutg[59] = getCodon("aug");

        cutg[60] = getCodon("ugg");

        cutg[61] = getCodon("uaa");
        cutg[62] = getCodon("uag");
        cutg[63] = getCodon("uga");
    }

    private static String stringifyCodon(BasisSymbol codon)
        throws IllegalSymbolException, BioException
    {
        // get the component symbols
        List codonList = codon.getSymbols();

        // get a tokenizer
        SymbolTokenization toke = RNATools.getRNA().getTokenization("token");

        String tokenizedCodon = toke.tokenizeSymbol((Symbol) codonList.get(0))
            + toke.tokenizeSymbol((Symbol) codonList.get(1))
            + toke.tokenizeSymbol((Symbol) codonList.get(2));

        return tokenizedCodon;
    }

    /**
     * writes out a CodonPref object in XML form
     */
    static void dumpToXML(CodonPref codonPref, XMLWriter xw, boolean writeWrapper)
        throws NullPointerException, IOException, IllegalSymbolException, BioException
    {
        // validate both objects first
        if ((codonPref == null) || (xw == null))
            throw new NullPointerException();

        // get the CodonPref Distribution
        Distribution codonDist = codonPref.getFrequency();

        // start <CodonPrefs>
        if (writeWrapper) xw.openTag("CodonPrefs");

        xw.openTag("CodonPref");
        xw.attribute("id", codonPref.getName());
        xw.attribute("geneticCodeId", codonPref.getGeneticCodeName());

        // loop over all codons, writing out the stats
        for (Iterator codonI = RNATools.getCodonAlphabet().iterator(); codonI.hasNext(); ) {
            BasisSymbol codon = (BasisSymbol) codonI.next();

            xw.openTag("frequency");

            // convert codon to a three letter string
            xw.attribute("codon", stringifyCodon(codon));
            xw.attribute("value", Double.toString(codonDist.getWeight(codon)));

            xw.closeTag("frequency");
        }

        xw.closeTag("CodonPref");

        if (writeWrapper) xw.closeTag("CodonPrefs");
    }

    /**
     * reads in records in CUD format
     */
    private static void readFromCUD(BufferedReader rdr, CodonPrefFilter filter)
    {
        try {
            String currLine;
            while ((currLine = rdr.readLine()) != null) {

                // process comment line
                StringTokenizer toke = new StringTokenizer(currLine, ":");
                if (toke.hasMoreTokens()) {
                    // get id string
                    String id = (toke.nextToken()).trim();

                    // read the codon count
                    currLine = rdr.readLine();
                    if (currLine == null) break;

                    // do we even want to process this record?
                    if (filter.isRequired(id)) {
                        toke = new StringTokenizer(currLine);

                        int idx = 0;
                        IndexedCount count = new IndexedCount(RNATools.getCodonAlphabet());
                        while (toke.hasMoreTokens()) {
                            // check that I haven't read too many values!
                            if (idx > 63) continue;
                            count.increaseCount(cutg[idx], Double.parseDouble(toke.nextToken()));
                            idx++;
                        }

                        if (idx != 64) continue;

                        // ok, I now have the counts and the name, let's stash it
                        Distribution codonDist = DistributionTools.countToDistribution(count);

                        CodonPref codonPref = new SimpleCodonPref("UNIVERSAL", codonDist, id);
                        filter.put(codonPref);
                    }
                }
            }
        }
        catch (IOException ioe) {}
        catch (IllegalSymbolException ise) {}
        catch (IllegalAlphabetException iae) {}
        catch (ChangeVetoException cve) {}
        catch (BioException be) {}
    }
}

