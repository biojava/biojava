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

package org.biojava.bio.seq.io;

import org.biojava.bio.BioException;
import org.biojava.bio.SimpleAnnotation;
import org.biojava.bio.seq.Feature;
import org.biojava.bio.seq.StrandedFeature;

/**
 * Simple parser for feature tables. This is shared between the EMBL
 * and GENBANK format readers.
 *
 * @author Thomas Down
 * @author Matthew Pocock
 * @author Greg Cox
 * @author Keith James
 * @deprecated Use org.biojavax.bio.seq.io framework instead
 */

/*
 * Greg Cox: Changed private fields and methods to protected so that
 *           SwissProtFeatureTableParser could subclass and snag the
 *           implementation.
 *
 * Thomas Down: Post 1.1, finally got round to refactoring this to be
 *              a `nice' player in the newio world.  Needless to say,
 *              this simplified things quite a bit.
 *
 * Keith James: Added support for reading fuzzy i.e. (123.567)
 *              locations in addition to unbounded i.e. <123..567
 *              locations.
 */

public class FeatureTableParser {
    private final static int   WITHOUT = 0;
    private final static int    WITHIN = 1;
    private final static int  LOCATION = 2;
    private final static int ATTRIBUTE = 3;

    private int featureStatus = WITHOUT;
    private StringBuffer featureBuf;
    private Feature.Template featureTemplate;

    private String                 featureSource;
    private SeqIOListener          listener;
    private EmblLikeLocationParser locParser;
    private String                 seqID;

    FeatureTableParser(SeqIOListener listener, String source) {
        this.listener      = listener;
        this.featureSource = source;
        //this.seqID = seqID;

        featureBuf = new StringBuffer();
        locParser  = new EmblLikeLocationParser(seqID);
    }

    public void setSeqID(String seqID) {
      this.seqID = seqID;
    }

    //
    // Interface which the processors use to call us
    //

    public void startFeature(String type) throws BioException {
        featureStatus = LOCATION;
        featureBuf.setLength(0);

        if (this.featureSource.equals("RefSeq:Protein")) {
            featureTemplate= new Feature.Template();
        }
        else {
            featureTemplate = new StrandedFeature.Template();
        }
        featureTemplate.type = type;
        featureTemplate.source = featureSource;
        featureTemplate.annotation = new SimpleAnnotation();
    }

    public void featureData(String line) throws BioException {
        switch (featureStatus) {
            case LOCATION:
                featureBuf.append(line);
                if (countChar(featureBuf, '(') == countChar(featureBuf, ')')) {
                    featureTemplate = locParser.parseLocation(featureBuf.substring(0), featureTemplate);
                    listener.startFeature(featureTemplate);
                    featureStatus = WITHIN;
                }
                break;

            case WITHIN:
                if (line.charAt(0) == '/') {
                    // System.out.println("got '/', quotes = " + countChar(line, '"'));
                    // attribute either is unquoted and on one line or
                    // is quoted, and must start & end with a quote
                    //
                    // we assume that no attributes have embedded quotes
                    int eq = line.indexOf("=");
                    if (line.charAt(eq + 1) != '"' ||
                        line.charAt(line.length() - 1) == '"'
                    ) {
                        processAttribute(line);
                    } else {
                        featureBuf.setLength(0);
                        featureBuf.append(line);
                        featureStatus = ATTRIBUTE;
                    }
                } else {
                    throw new BioException("Invalid line in feature body: " + line);
                }
                break;

            case ATTRIBUTE:
                // If the attribute contains whitespace it probably
                // consists of whitespace-delimited words. Therefore a
                // space should be inserted at EOL otherwise words will
                // get fused (unless there is a space already there)
                if (((featureBuf.toString().indexOf(" ") >= 0) ||
                     (line.toString().indexOf(" ") >= 0)) &&
                    featureBuf.toString().charAt(featureBuf.length()-1) != ' '){
                    featureBuf.append(" ");
                }
                featureBuf.append(line);
                

                int eq = featureBuf.toString().indexOf("=");
                if (featureBuf.charAt(eq + 1) != '"' ||
                    featureBuf.charAt(featureBuf.length() - 1) == '"'
                ) {
                    processAttribute(featureBuf.substring(0));
                    featureStatus = WITHIN;
                }
                break;
        }
    }

    public void endFeature()
        throws BioException {
        listener.endFeature();
        featureStatus = WITHOUT;
    }

    public boolean inFeature() {
        return (featureStatus != WITHOUT);
    }

    /**
     * Process the a string corresponding to a feature-table
     * attribute, and fire it off to our listener.
     */
    private void processAttribute(String attr) throws BioException {
        // System.err.println(attr);
        int eqPos = attr.indexOf('=');
        if (eqPos == -1) {
            listener.addFeatureProperty(attr.substring(1), Boolean.TRUE);
        } else {
            String tag = attr.substring(1, eqPos);
            eqPos++;

            if (attr.charAt(eqPos) == '"')
                ++eqPos;
            int max = attr.length();

            if (attr.charAt(max - 1) == '"')
                --max;
            String val = attr.substring(eqPos, max);

            if (val.indexOf('"') >= 0) {
                StringBuffer sb = new StringBuffer();
                boolean escape = false;
                for (int i = 0; i < val.length(); ++i) {
                    char c = val.charAt(i);
                    if (c == '"') {
                        if (escape)
                            sb.append(c);
                        escape = !escape;
                    } else {
                        sb.append(c);
                        escape = false;
                    }
                }
                val = sb.substring(0);
            }
            listener.addFeatureProperty(tag, val);
        }
    }

    private int countChar(StringBuffer s, char c) {
        int cnt = 0;
        int length = s.length();
        for (int i = 0; i < length; ++i)
            if (s.charAt(i) == c)
                ++cnt;
        return cnt;
    }
}
