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

/*
 *    AAindexStreamReader.java
 */
package org.biojava.bio.proteomics.aaindex;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.Reader;
import java.util.Map;
import java.util.NoSuchElementException;

import org.biojava.bio.BioException;
import org.biojava.bio.seq.io.SymbolTokenization;
import org.biojava.bio.symbol.Symbol;
import org.biojava.bio.symbol.SymbolPropertyTable;

/**
 * Iterator over {@link org.biojava.bio.proteomics.aaindex.AAindex} objects that 
 * are stored in a stream in the AAindex1 file format. The format
 * of such an Amino Acid Index Database file is described in the
 * <a href="http://www.genome.ad.jp/dbget-bin/show_man?aaindex">AAindex manual
 * </a>. The {@link #nextTable()} method returns objects of type
 * {@link org.biojava.bio.proteomics.aaindex.AAindex}. See this class also for
 * further informations. To hold an AAindex1 file in memory for random access
 * use the {@link org.biojava.bio.proteomics.aaindex.SimpleSymbolPropertyTableDB}
 * class:
 * <pre>
 * SimpleSymbolPropertyTableDB db = new SimpleSymbolPropertyTableDB(
 *         new AAindexStreamReader(new FileReader("aaindex1")));
 * AAindex hydrophobicity = (AAindex) db.table("CIDH920105");
 * SymbolList symbols = ProteinTools.createProtein(
 *     "ARNDCEQGHILKMFPSTWYV");
 * double hp = 0.0;
 * for (int i = 1; i <= symbols.length(); i++) {
 *     hp += hydrophobicity.getDoubleValue(symbols.symbolAt(i));
 * }
 * System.out.println("Average hydrophobicity: " + Double.toString(
 *         hp / symbols.length()));
 * </pre>
 * <p><b>References:</b></p>
 * <p><a href="http://www.genome.ad.jp/dbget/aaindex.html">AAindex web 
 * site</a>.</p>
 * <p>Kawashima, S. and Kanehisa, M.; AAindex: amino acid index database. 
 * Nucleic Acids Res. 28, 374 (2000).</p>
 * <p>Tomii, K. and Kanehisa, M.;  Analysis of amino acid indices and mutation 
 * matrices for sequence comparison and structure prediction of proteins. 
 * Protein Eng. 9, 27-36 (1996).</p>
 * <p>Nakai, K., Kidera, A., and Kanehisa, M.;  Cluster analysis of amino acid 
 * indices for prediction of protein structure and function.  
 * Protein Eng. 2, 93-100 (1988)</p>
 * @author <a href="mailto:Martin.Szugat@GMX.net">Martin Szugat</a>
 * @version $Revision$
 */
public class AAindexStreamReader implements SymbolPropertyTableIterator {
    
//    public static final void main(String[] args) throws NullPointerException, 
//    FileNotFoundException, BioException, IOException {
//        SimpleSymbolPropertyTableDB db = new SimpleSymbolPropertyTableDB(
//                new AAindexStreamReader(new FileReader("aaindex1")));
//        AAindex hydrophobicity = (AAindex) db.table("CIDH920105");
//        SymbolList symbols = ProteinTools.createProtein(
//            "ARNDCEQGHILKMFPSTWYV");
//        double hp = 0.0;
//        for (int i = 1; i <= symbols.length(); i++) {
//            hp += hydrophobicity.getDoubleValue(symbols.symbolAt(i));
//        }
//        System.out.println("Average hydrophobicity: " + Double.toString(
//                hp / symbols.length()));
//    }
//        
    /* PRIVATE CONSTANTS */

    /**
     * Name of the tokenizer.
     */
    private static final String TOKENIZER = "token";

    /* STATIC FIELDS */

    /**
     * List of amino acid symbols.
     */
    private static Symbol[] aa = null;

    /* STATIC CONSTRUCTOR */

    static {
        try {
            SymbolTokenization tokenizer = 
                AAindex.PROTEIN_ALPHABET.getTokenization(TOKENIZER);
            aa = new Symbol[] {tokenizer.parseToken("A"),
                    tokenizer.parseToken("R"), tokenizer.parseToken("N"),
                    tokenizer.parseToken("D"), tokenizer.parseToken("C"),
                    tokenizer.parseToken("Q"), tokenizer.parseToken("E"),
                    tokenizer.parseToken("G"), tokenizer.parseToken("H"),
                    tokenizer.parseToken("I"), tokenizer.parseToken("L"),
                    tokenizer.parseToken("K"), tokenizer.parseToken("M"),
                    tokenizer.parseToken("F"), tokenizer.parseToken("P"),
                    tokenizer.parseToken("S"), tokenizer.parseToken("T"),
                    tokenizer.parseToken("W"), tokenizer.parseToken("Y"),
                    tokenizer.parseToken("V"), };
        } catch (BioException e) {
            e.printStackTrace();
        } catch (IndexOutOfBoundsException e) {
            e.printStackTrace();
        }
    };
    
    /* PRIVATE FIELDS */

    /**
     * The internal reader.
     */
    private BufferedReader reader = null;

    /**
     * The current read line.
     */
    private String line = null;

    /**
     * The key char of the current read section. 
     */
    private char keyChar;

    /**
     * The value of the current read section.
     */
    private String stringValue;
    
    /* PUBLIC CONSTRUCTORS */

    /**
     * Initializes the iterator.
     * @param reader reader over a stream in the AAindex file format.
     * @throws IOException if the stream could not be read.
     * @throws NullPointerException if <code>reader</code> is <code>null</code>.
     */
    public AAindexStreamReader(Reader reader) throws IOException, 
    NullPointerException {
        this(new BufferedReader(reader));
    }

    /**
     * Initializes the iterator.
     * @param reader buffered reader over a stream in the AAindex file format.
     * @throws IOException if the stream could not be read.
     * @throws NullPointerException if <code>reader</code> is <code>null</code>.
     */
    public AAindexStreamReader(BufferedReader reader) throws IOException,
    NullPointerException {
        if (reader == null) {
            throw new NullPointerException("reader is null.");
        }
        this.reader = reader;
        line = reader.readLine();
    }
    
    /* PUBLIC METHODS */
    
    /**
     * Checks if the end of the file or stream is reached.
     * @return <code>true</code> if the end of the file is reached,
     * <code>false</code> otherwise. 
     */
    public boolean eof() {
        if (line == null) {
            return true;
        } else {
            while (line != null && line.length() == 0) {
                try {
                    line = reader.readLine();
                } catch (IOException e) {
                    return true;
                }
            }
            return (line == null);
        }
    }

    /**
     * Reads a AAindex section.
     * @throws BioException if the section could not be read.
     */
    private void readSection() throws BioException {

        keyChar = line.charAt(0);

        StringBuffer stringBuffer = new StringBuffer();

        do {
            if (line.length() > 2) {
                stringBuffer.append(line.substring(2));
                if (!line.endsWith(" ")) {
                    stringBuffer.append(" ");
                }
            }
            try {
                line = reader.readLine();
            } catch (IOException e) {
                throw new BioException(e);
            }
        } while (!eof() && line.charAt(0) == ' ');

        stringValue = stringBuffer.toString();
    }

    /* INTERFACE SymbolPropertyTableIterator */

    /**
     * {@inheritDoc}
     */
    public boolean hasNext() {
        return (!eof());
    }

    /**
     * {@inheritDoc}
     */
    public SymbolPropertyTable nextTable() throws BioException {

        if (eof()) {
            throw new NoSuchElementException();
        }

        readSection();

        if (keyChar != 'H') {
            throw new BioException("Expected 'H' but found: '" + keyChar
                    + "'.");
        }
        AAindex aaIndex = new AAindex(stringValue.trim());

        readSections: while (!eof()) {

            readSection();

            switch (keyChar) {
            case 'D':
                aaIndex.setDescription(stringValue);
                break;
            case 'R':
                aaIndex.setLITDBEntryNumbers(stringValue.split("\\s+"));
                break;
            case 'A':
                aaIndex.setArticleAuthors(stringValue);
                break;
            case 'T':
                aaIndex.setArticleTitle(stringValue);
                break;
            case 'J':
                aaIndex.setJournalReference(stringValue);
                break;
            case 'C':
                String[] keyValuePairs = stringValue.split("\\s+");
                Map similarEntries = aaIndex.similarEntries();
                for (int i = 0; i < keyValuePairs.length - 1; i += 2) {
                    similarEntries.put(keyValuePairs[i], Double
                            .valueOf(keyValuePairs[i + 1]));
                }
                break;
            case 'I':
                String[] headersAndIndices = stringValue.split("\\s+");
                for (int i = 0; i < 20; i++) {
                    try {
                        aaIndex.setDoubleProperty(aa[i],
                                headersAndIndices[11 + i]);
                    } catch (NumberFormatException e) {
                        aaIndex.setDoubleProperty(aa[i],
                                        "NaN");
                    }
                }
                break;
            case '*':
                aaIndex.setComment(stringValue);
                break;
            case '/':
                break readSections;
            default:
                throw new BioException("Invalid key char found: " + keyChar
                        + "'.");
            }
        }
        
        return aaIndex;
    }
}
