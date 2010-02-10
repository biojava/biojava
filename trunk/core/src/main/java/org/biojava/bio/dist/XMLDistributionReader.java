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
 */


package org.biojava.bio.dist;

import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;

import javax.xml.parsers.ParserConfigurationException;
import javax.xml.parsers.SAXParser;
import javax.xml.parsers.SAXParserFactory;

import org.biojava.bio.BioException;
import org.biojava.bio.seq.io.SymbolTokenization;
import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.AlphabetManager;
import org.biojava.bio.symbol.AtomicSymbol;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.Symbol;
import org.biojava.utils.ChangeVetoException;
import org.xml.sax.Attributes;
import org.xml.sax.InputSource;
import org.xml.sax.SAXException;
import org.xml.sax.helpers.DefaultHandler;


/**
 * A SAX parser that reads an XML representation of a
 * Distribution from a file and recreates it as a Distribution Object.
 * Handles OrderNDistributions and Simple Distributions but
 * ensure the OrderNDistributions being read in was made using
 * conditioning and conditioned Alphabets.
 *
 * @author Russell Smithies
 * @author Mark Schreiber
 * @author Matthew Pocock
 * @since 1.3
 */
public class XMLDistributionReader extends DefaultHandler {
  //fixme: the handler and the user API seem to be muddled up - can the
  // DefaultHandler impl be factored out into a static and hopefully private
  // class?
    private Alphabet alpha = null;
    private Distribution dist = null;
    private DistributionFactory fact = null;
    //private OrderNDistributionFactory ondFact = null;
    private SymbolTokenization nameParser = null;
    private Symbol sym = null;

    private Alphabet conditioningAlpha = null;
    private AtomicSymbol conditioningSymbol = null;
    private SymbolTokenization conditioningTok = null;
    private Alphabet conditionedAlpha = null;
    private AtomicSymbol conditionedSymbol = null;
    private SymbolTokenization conditionedTok = null;

    private Distribution getDist() {
        return dist;
    } //end getDist

    /**
     * Reads an XML representation of a Distribution from a file.
     *
     * @param is input in XML format
     * @return dist the Distribution created.
     * @throws IOException if an error occurs during reading.
     * @throws SAXException if the XML is not as expected.
     */
    public Distribution parseXML(InputStream is) throws IOException, SAXException{
        //org.xml.sax.XMLReader parser = new org.apache.xerces.parsers.SAXParser();
        SAXParserFactory fact = SAXParserFactory.newInstance();
        fact.setNamespaceAware(true);
        try {
          SAXParser parser = fact.newSAXParser();

//          parser.setContentHandler(this);
//
//          parser.setErrorHandler(this);

          InputSource xml = null;

          xml = new InputSource(new InputStreamReader(is));

          parser.parse(xml, this);
        }
        catch (ParserConfigurationException ex) {
         //really shouldn't happen
         throw new SAXException("Cannot make SAXParser",ex);
        }

        return this.getDist();
    }

    /**
     * Required by SAXParser to be public.
     *
     * <p>
     * It is not reccomended that you use this method directly. Use ParseXML
     * instead.
     * </p>
     */
    public void startElement(String nameSpaceURI, String localName, String rawName, Attributes attributes)
                      throws SAXException{
        if (localName.equals("Distribution") || localName.equals("OrderNDistribution")) {
            processDistElement(attributes);
        } else if (localName.equals("alphabet")) {
            processAlphabetElement(attributes);
        } else if (localName.equals("conditioning_symbol")) {
            processConditioningSymbol(attributes);
        } else if (localName.equals("weight")) {
            processWeightElement(attributes);
        }
    } //end startElement

    private void processConditioningSymbol(Attributes attr)
                                    throws SAXException {
        String name = attr.getValue("name");
        try {
          conditioningSymbol = (AtomicSymbol)conditioningTok.parseToken(name);
        }
        catch (IllegalSymbolException ex) {
          throw new SAXException(ex);
        }

    }

    private void processWeightElement(Attributes attr)
                               throws SAXException{
        double weight = 0.0;

        try {
            //get the weight of the symbol
            weight = Double.parseDouble(attr.getValue("prob"));
        } catch (NumberFormatException ex) {
            //catches the "NAN" string
            weight = 0.0;
        }

        //add counts if SIMPLE DISTRIBUTION
        if ((dist instanceof OrderNDistribution) == false) {
            try {
                //initialize the tokenizer
                nameParser = dist.getAlphabet().getTokenization("name");
            } catch (BioException ex) {
                throw new SAXException("Couldn't get tokenization for "
                                       +dist.getAlphabet().getName(), ex);
            }

            try {
                //get the symbol name
                sym = nameParser.parseToken(attr.getValue("sym"));


                //add count to dist
                dist.setWeight(sym, weight);
            } catch (IllegalSymbolException ex) {
                throw new SAXException("Illegal symbol found", ex);
            } catch (ChangeVetoException ex) {
                throw new SAXException("Distribution has been locked, possible synchronization problem !?",ex);
            }

            //add countf if ORDER N DISTRIBUTION
        } else if (dist instanceof OrderNDistribution) {


            //get the weight for symbol
            try {
                //get the weight of the symbol
                weight = Double.parseDouble(attr.getValue("prob"));
            } catch (NumberFormatException ex) {
                weight = 0.0;
            }

            //rebuild the symbol from the conditioning and conditioned symbol
            String name = attr.getValue("sym");
            try {
              conditionedSymbol = (AtomicSymbol)conditionedTok.parseToken(name);
            }
            catch (IllegalSymbolException ex) {
              throw new SAXException(ex);
            }
            List l = new ArrayList();
            l.add(conditioningSymbol);
            l.add(conditionedSymbol);


            try {
                sym = alpha.getSymbol(l);
                //set weights on distribution
                dist.setWeight(sym, weight);
            } catch (IllegalSymbolException ex) {
                throw new SAXException("Illegal symbol found", ex);
            } catch (ChangeVetoException ex) {
                throw new SAXException("Distribution has been locked, possible synchronization problem !?",ex);
            }
        }
    } //end processWeight

    private void processAlphabetElement(Attributes attr)
                                 throws SAXException {
        String alphaName = attr.getValue("name");


        //get Alphabet
        alpha = AlphabetManager.alphabetForName(alphaName);

        //make the Distribution
        try {
            dist = fact.createDistribution(alpha);
            if(dist instanceof OrderNDistribution){
              conditionedAlpha = ((OrderNDistribution)dist).getConditionedAlphabet();
              conditionedTok = conditionedAlpha.getTokenization("name");

              conditioningAlpha = ((OrderNDistribution)dist).getConditioningAlphabet();
              conditioningTok = conditioningAlpha.getTokenization("name");
            }

        } catch (IllegalAlphabetException ex) {
            throw new SAXException(ex);
        } catch (BioException ex) {
            throw new SAXException(ex);
        }
    } //end processAlphabetElement

    private void processDistElement(Attributes attr) throws SAXException {
        if (attr.getValue("type").equals("Distribution")) {
            fact = DistributionFactory.DEFAULT;
        } else if (attr.getValue("type").equals("OrderNDistribution")) {
            fact = OrderNDistributionFactory.DEFAULT;
        } else {
            throw new SAXException("Element must be a distribution");
        }
    } //end processDistElement

}