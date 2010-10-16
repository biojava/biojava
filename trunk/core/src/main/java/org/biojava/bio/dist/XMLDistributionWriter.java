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

import java.io.BufferedWriter;
import java.io.IOException;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.util.Iterator;

import org.biojava.bio.BioError;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.Symbol;


/**
 * Writes an OrderNDistribution or simple Distribution to an XML file.
 *
 * @author Russell Smithies
 * @author Mark Schreiber
 * @since 1.3
 */
public class XMLDistributionWriter{

  BufferedWriter out = null;


    private void writeOND2XML(OrderNDistribution ond, OutputStream os) throws IOException{
        Distribution conditionedDist = null;
        FiniteAlphabet conditionedAlpha = null;
        out = new BufferedWriter(new OutputStreamWriter(os));


        out.write("<?xml version=\"1.0\" ?>");
        out.write("<Distribution type=\"OrderNDistribution\">");
        out.write("<alphabet name=\"" + ond.getAlphabet().getName() +
                  "\"/>");

        for (Iterator i = ((FiniteAlphabet) ond.getConditioningAlphabet()).iterator();
             i.hasNext();) {
          Symbol sym = (Symbol) i.next();
          out.write("<conditioning_symbol name=\"" + sym.getName() +
                    "\">");

          try {
            conditionedDist = ond.getDistribution(sym);
          } catch (IllegalSymbolException ex) {
            throw new BioError("Distribution has been built with Illegal Symbols !?", ex);
          }

          conditionedAlpha = (FiniteAlphabet) conditionedDist.getAlphabet();

          for (Iterator it = conditionedAlpha.iterator(); it.hasNext();) {
            Symbol condSym = (Symbol) it.next();
            double weight = 0.0;

            try {
              weight = conditionedDist.getWeight(condSym);
            } catch (IllegalSymbolException ex) {
              throw new BioError("Distribution has been built with Illegal Symbols !?", ex);
            }

            out.write("<weight sym=\"" + condSym.getName() +
                      "\" prob=\"" + weight + "\"/>");
          }

          out.write("</conditioning_symbol>");
        }

        out.write("</Distribution>");
        out.flush();
    } //end writeXML


    private void writeDist2XML(Distribution d, OutputStream os) throws IOException{
         out = new BufferedWriter(new OutputStreamWriter(os));


         out.write("<?xml version=\"1.0\" ?>");
         out.write("<Distribution type=\"Distribution\">");
         out.write("<alphabet name=\"" + d.getAlphabet().getName() +
                   "\"/>");

         for (Iterator i = ((FiniteAlphabet) d.getAlphabet()).iterator();
              i.hasNext();) {
           Symbol sym = (Symbol) i.next();
           double weight = 0.0;

           try {
             weight = d.getWeight(sym);
           } catch (IllegalSymbolException ex) {
             throw new BioError("Distribution has been built with Illegal Symbols !?", ex);
           }

           out.write("<weight sym=\"" + sym.getName() +
                     "\" prob=\"" + weight + "\"/>");
         }
         out.write("</Distribution>");
         out.flush();

    } //end writeXML

    /**
     * Writes an OrderNDistribution or simple Distribution to an XML file.
     * @param d  Distribution to write.
     * @param os OutputStream to write to.
     * @throws IOException if an error occurs during writing
     */
    public void writeDistribution(Distribution d, OutputStream os) throws IOException{
        if (d instanceof OrderNDistribution) {
            writeOND2XML((OrderNDistribution) d, os);
        } else {
            writeDist2XML(d, os);
        }

    }
}