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

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.StringTokenizer;

import org.biojava.bio.Annotation;
import org.biojava.bio.BioException;
import org.biojava.utils.ChangeVetoException;
import org.biojava.utils.CommitFailure;
import org.biojava.utils.SmallMap;
import org.biojava.utils.io.RAF;

/**
 * <code>BioStore</code>s represent directory and file structures
 * which index flat files according to the OBDA specification. The
 * preferred method of constructing new instances is to use
 * <code>BioStoreFactory</code>.
 *
 * @author Matthew Pocock
 * @author Keith James
 */
public class BioStore implements IndexStore {

    /**
     * <code>STRING_CASE_SENSITIVE_ORDER</code> compares two
     * <code>Object</code>s, which must both be <code>String</code>s,
     * lexicographically using <code>compareTo</code>. The comparison
     * is carried out 'a' to 'b'.
     */
    static Comparator STRING_CASE_SENSITIVE_ORDER = new Comparator() {
            public int compare(Object a, Object b) {
                return ((Comparable) a).compareTo(b);
            }
        };
    
    private ConfigFile metaData;
    private File location;
    private String primaryKey;
    private Map idToList;
    private RAF[] fileIDToRAF;
    private SearchableList primaryList;
    private int fileCount;

    /**
     * Creates a new <code>BioStore</code> flatfile index at the
     * specified location with the specified caching behaviour.
     *
     * @param location a <code>File</code> indicating the index
     * directory.
     * @param cache a <code>boolean</code> indicating whether the
     * implementation should cache its state.
     *
     * @exception IOException if an error occurs.
     * @exception BioException if an error occurs.
     */
    public BioStore(File location, boolean cache)
        throws IOException, BioException {
        this(location, cache, false);
    }
    
    BioStore(File location, boolean cache, boolean mutable)
        throws IOException, BioException {
        this.location = location;

        File configFile = BioStoreFactory.makeConfigFile(location);
        if (!configFile.exists()) {
            throw new BioException("Config file does not exist: "
                                   + configFile);
        }
        metaData = new ConfigFile(BioStoreFactory.makeConfigFile(location));
        idToList = new SmallMap();

        primaryKey = (String) metaData.getProperty(BioStoreFactory.PRIMARY_KEY_NAME);
        String keyList = (String) metaData.getProperty(BioStoreFactory.KEYS);

        File plFile = BioStoreFactory.makePrimaryKeyFile(location, primaryKey);
        if (cache) {
            primaryList = new CacheList(new PrimaryIDList(plFile, this, mutable));
        } else {
            primaryList = new PrimaryIDList(plFile, this, mutable);
        }

        StringTokenizer sTok = new StringTokenizer(keyList, "\t");
        while (sTok.hasMoreTokens()) {
            String k = sTok.nextToken();

            File file = BioStoreFactory.makeSecondaryFile(location, k);
            if (cache) {
                idToList.put(k, new CacheList(new SecondaryFileAsList(file, mutable)));
            } else {
                idToList.put(k, new SecondaryFileAsList(file, mutable));
            }
        }

        readFileIDs();
    }
    
    /**
     * The name of this store or null if the name has not been set.
     */
    public String getName() {
      if (metaData.containsProperty(BioStoreFactory.STORE_NAME)) {
        return (String) metaData.getProperty(BioStoreFactory.STORE_NAME);
      } else {
        return null;
      }
    }
    
    /**
     * <code>getLocation</code> returns the directory where the index
     * is located.
     *
     * @return a <code>File</code>.
     */
    public File getLocation() {
      return location;
    }

    private void readFileIDs()
        throws
            IOException,
            BioException
    {
        fileIDToRAF = new RAF[5];
        fileCount = 0;

        for (Iterator i = metaData.keys().iterator(); i.hasNext(); ) {
            String key = (String) i.next();
            if (key.startsWith("fileid_")) {
                int indx = Integer.parseInt(key.substring("fileid_".length()));
                String fileLine = (String) metaData.getProperty(key);
                int tab = fileLine.indexOf("\t");
                File file = new File(fileLine.substring(0, tab));
                RAF raf = new RAF(file, "r");
                long length = Long.parseLong(fileLine.substring(tab+1));

                if (file.length() != length) {
                    throw new BioException("File changed length: " + file);
                }
                
                if (indx >= fileCount) {
                    // beyond end
                    if (indx >= fileIDToRAF.length) {
                        // beyond array end
                        RAF[] tmpr = new RAF[indx+1];
                        System.arraycopy(fileIDToRAF, 0, tmpr, 0, fileIDToRAF.length);
                        fileIDToRAF = tmpr;
                    }

                    fileCount = indx;
                }
                fileIDToRAF[indx] = raf;
            }
        }
    }

    private void writeFileIDs()
        throws BioException, IOException, ChangeVetoException {
        for (int i = 0; i < fileCount; i++) {
            RAF file = fileIDToRAF[i];
            long length = file.length();
            String prop = "fileid_" + i;
            String val = file.getFile().toString() + "\t" + length;
            metaData.setProperty(prop, val);
        }
    }

    RAF getFileForID(int fileId) {
        return fileIDToRAF[fileId];
    }

    int getIDForFile(RAF file)
        throws IOException {
        // scan list
        for (int i = 0; i < fileCount; i++) {
            if (file.equals(fileIDToRAF[i])) {
                return i;
            }
        }

        // extend fileIDToFile array
        if (fileCount >= fileIDToRAF.length) {
            RAF[] tmpr = new RAF[fileIDToRAF.length + 4]; // 4 is magic number
            System.arraycopy(fileIDToRAF, 0, tmpr, 0, fileCount);
            fileIDToRAF = tmpr;
        }

        // add the unseen file to the list
        fileIDToRAF[fileCount] = file;
        return fileCount++;
    }

    public Annotation getMetaData() {
        return metaData;
    }

    public Record get(String id) {
        return (Record) primaryList.search(id);
    }

    public List get(String id, String namespace)
        throws BioException {
        List hits = new ArrayList();
        if (namespace.equals(primaryKey)) {
            hits.add(primaryList.search(id));
        } else {
            SecondaryFileAsList secList = (SecondaryFileAsList) idToList.get(namespace);
            List kpList = secList.searchAll(id);
            for (Iterator i = kpList.iterator(); i.hasNext(); ) {
                KeyPair keyPair = (KeyPair) i.next();
                hits.add(primaryList.search(keyPair.getSecondary()));
            }
        }

        return hits;
    }

    public void writeRecord(RAF file,
                            long offset,
                            int length,
                            String id,
                            Map secIDs) {
        primaryList.add(new Record.Impl(id, file, offset, length));
        if (!secIDs.isEmpty()) {
            for (Iterator mei = secIDs.entrySet().iterator(); mei.hasNext(); ) {
                Map.Entry me = (Map.Entry) mei.next();
                String sid = (String) me.getKey();
                List sfl = (List) idToList.get(sid);
                Collection svals = (Collection) me.getValue();
                for (Iterator i = svals.iterator(); i.hasNext(); ) {
                    String sval = (String) i.next();
                    sfl.add(new KeyPair.Impl(sval, id));
                }
            }
        }
    }

    /**
     * <code>getRecordList</code> returns all the <code>Record</code>s
     * in the index.
     *
     * @return a <code>List</code> of <code>Record</code>s.
     */
    public List getRecordList() {
        return primaryList;
    }

    /**
     * <code>commit</code> writes an index to disk.
     *
     * @exception CommitFailure if an error occurs.
     */
    public void commit()
        throws CommitFailure {
        Collections.sort(primaryList, primaryList.getComparator());
        primaryList.commit();
        for (Iterator i = idToList.values().iterator(); i.hasNext(); ) {
            SearchableList fal = (SearchableList) i.next();
            Collections.sort(fal, fal.getComparator());
            ((SearchableList) fal).commit();
        }

        try {
            writeFileIDs();
        } catch (ChangeVetoException cve) {
            throw new CommitFailure(cve);
        } catch (IOException ioe) {
            throw new CommitFailure(ioe);
        } catch (BioException be) {
          throw new CommitFailure(be);
        }

        metaData.commit();
    }
}
