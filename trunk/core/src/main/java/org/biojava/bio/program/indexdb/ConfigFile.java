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

package org.biojava.bio.program.indexdb;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

import org.biojava.bio.AbstractAnnotation;
import org.biojava.utils.AssertionFailure;
import org.biojava.utils.CommitFailure;
import org.biojava.utils.Commitable;

/**
 * <code>ConfigFile</code> implements reading, updating and writing
 * of OBDA tab-delimited index files.
 *
 * @author Matthew Pocock
 * @author Keith James
 */
class ConfigFile
    extends AbstractAnnotation
    implements Commitable
{
    private File file;
    private Map map;

    /**
     * Creates a new <code>ConfigFile</code> using the specified
     * <code>File</code>, which may or may not already exist.
     *
     * @param file a <code>File</code>.
     *
     * @exception IOException if an error occurs.
     */
    public ConfigFile(File file)
        throws IOException {
        this.file = file;
        map = new HashMap();
        if(file.exists()) {
            parseFile();
        }
    }

    public void commit()
        throws CommitFailure {
        try {
            writeFile();
        } catch (IOException e) {
            throw new CommitFailure("Couldn't commit", e);
        }
    }

    public void rollback() {
        try {
            parseFile();
        } catch (IOException e) {
            throw new AssertionFailure("Couldn't roll back: your data may be invalid", e);
        }
    }

    private void parseFile()
        throws IOException {
        BufferedReader reader = new BufferedReader(new FileReader(file));

        for(String line = reader.readLine();
            line != null;
            line = reader.readLine()) {
            int tab = line.indexOf("\t");
            String key = line.substring(0, tab).trim();
            String value = line.substring(tab+1).trim();

            map.put(key, value);
        }
    }

    private void writeFile()
        throws IOException {
        PrintWriter writer = new PrintWriter(new FileWriter(file));

        writer.println("index\t" + map.get("index"));

        for(Iterator i = map.entrySet().iterator(); i.hasNext(); ) {
            Map.Entry me = (Map.Entry) i.next();
            if(!me.getKey().equals("index")) {
                writer.println(me.getKey() + "\t" + me.getValue());
            }
        }

        writer.flush();
    }

    protected Map getProperties() {
        return map;
    }

    protected boolean propertiesAllocated() {
        return true;
    }
}
