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

import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;

import org.biojava.bio.BioException;
import org.biojava.utils.xml.PrettyXMLWriter;


public interface CodonPrefFilter
{
    /**
     * indicates if the current CodonPref is to be accepted
     */
    public boolean isRequired(String name);

    /**
     * handles storage of a CodonPref object
     */
    public void put(CodonPref codonPref) throws BioException;


    public class AcceptAll implements CodonPrefFilter {
      List filteredCodonPrefs;

      public AcceptAll(){
        filteredCodonPrefs = new ArrayList();
      }

      public boolean isRequired(String name){
        return true;
      }

      public void put(CodonPref codonPref){
        filteredCodonPrefs.add(codonPref);
      }

      public List getCodonPrefs(){
        return filteredCodonPrefs;
      }
    }

    public class ByName implements CodonPrefFilter
    {
        String name;
        CodonPref filteredCodonPref;

        public ByName(String name)
        {
            this.name = name;
        }

        public boolean isRequired(String name)
        {
            return name.equals(this.name);
        }

        public void put(CodonPref codonPref)
            throws BioException
        {
            filteredCodonPref = codonPref;
        }

        public CodonPref getCodonPref()
        {
            return filteredCodonPref;
        }
    }

    public class EverythingToXML implements CodonPrefFilter
    {
        PrintWriter pw;
        PrettyXMLWriter xw;
        boolean writtenWrapper = false;

        public EverythingToXML(PrintWriter pw)
        {
            this.pw = pw;
            xw = new PrettyXMLWriter(pw);
        }

        public boolean isRequired(String name)
        {
            return true;
        }

        public void put(CodonPref codonPref)
            throws BioException
        {
            try {
            if (!writtenWrapper) {
                xw.openTag("CodonPrefs");
                writtenWrapper = true;
            }

            CodonPrefTools.dumpToXML(codonPref, xw, false);
            }
            catch (IOException ioe) {
                throw new BioException(ioe);
            }
        }

        public void close()
            throws IOException
        {
            xw.closeTag("CodonPrefs");
            pw.flush();
        }
    }
}

