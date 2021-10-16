package org.biojava.nbio.structure.chem;

import org.biojava.nbio.core.util.InputStreamProvider;
import org.biojava.nbio.structure.align.util.URLConnectionTools;
import org.biojava.nbio.structure.align.util.UserConfiguration;
import org.biojava.nbio.structure.io.LocalPDBDirectory;
import org.biojava.nbio.structure.io.cif.ChemCompConverter;
import org.rcsb.cif.ParsingException;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.net.URL;
import java.net.URLConnection;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.nio.file.StandardCopyOption;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.zip.GZIPOutputStream;

/**
 * This provider of chemical components can download and cache chemical component definition files from the RCSB PDB web
 * site. It is the default way to access these definitions. If this provider is called he first time, it will download
 * and install all chemical component definitions in a local directory. Once the definition files have been installed,
 * it has quick startup time and low memory requirements.
 *
 * An alternative provider, that keeps all definitions in memory is the {@link AllChemCompProvider}. Another provider,
 * that does not require any network access, but only can support a limited set of chemical component definitions, is
 * the {@link ReducedChemCompProvider}.
 *
 * @author Andreas Prlic
 */
public class DownloadChemCompProvider implements ChemCompProvider {
    private static final Logger logger = LoggerFactory.getLogger(DownloadChemCompProvider.class);

    private static final String NEWLINE = System.getProperty("line.separator");

    public static final String CHEM_COMP_CACHE_DIRECTORY = "chemcomp";
    public static final String DEFAULT_SERVER_URL = "https://files.rcsb.org/ligands/download/";
    public static final String DEFAULT_CHEMCOMP_PATHURL_TEMPLATE = "{ccd_id}.cif";

    /**
     * The base URL to which the full path specified via {@link #setChemCompPathUrlTemplate(String)} is appended.
     * It is assumed that it has a trailing slash.
     */
    public static String serverBaseUrl = DEFAULT_SERVER_URL;

    private static File path;

    private static String chemCompPathUrlTemplate = DEFAULT_CHEMCOMP_PATHURL_TEMPLATE;

    static final Pattern CCD_ID_TEMPLATE_REGEX = Pattern.compile("\\{ccd_id(?::(\\d+_\\d+|[-+]?\\d+))?}");


    // flags to make sure there is only one thread running that is loading the dictionary
    static AtomicBoolean loading = new AtomicBoolean(false);

    static final List<String> protectedIDs = new ArrayList<>();
    static {
        protectedIDs.add("CON");
        protectedIDs.add("PRN");
        protectedIDs.add("AUX");
        protectedIDs.add("NUL");
    }

    private static ChemCompProvider fallback = null; // Fallback provider if the download fails

    /**
     * by default we will download only some of the files. User has to request that all files should be downloaded...
     */
    boolean downloadAll = false;

    public DownloadChemCompProvider() {
        this(null);
    }

    public DownloadChemCompProvider(String cacheFilePath) {
        logger.debug("Initialising DownloadChemCompProvider");

        // note that path is static, so this is just to make sure that all non-static methods will have path initialised
        if (cacheFilePath != null) {
            path = new File(cacheFilePath);
        }
    }

    /**
     * Set the base URL for the location of all chemical component CIF files, to which the chemCompPathUrlTemplate
     * is appended, settable in {@link #setChemCompPathUrlTemplate(String)}. A trailing slash is appended
     * if not present.
     */
    public static void setServerBaseUrl(String serverBaseUrl) {
        if (!serverBaseUrl.endsWith("/")) {
            serverBaseUrl = serverBaseUrl + "/";
        }
        DownloadChemCompProvider.serverBaseUrl = serverBaseUrl;
    }

    /**
     * Set the path to append to the serverBaseUrl (settable in {@link #setServerBaseUrl(String)}).
     * The string can contain placeholders that will be expanded at runtime:
     * <li>"{ccd_id}" to be replaced by the chemical component identifier, in capitals</li>
     * <li>"{ccd_id:beginIndex-endIndex}" to be replaced by a substring of the chemical component identifier in capitals,
     * with indices following the same convention as {@link String#substring(int, int)} </li>
     * <li>"{ccd_id:index}" to be replaced by a substring of the chemical component identifier in capitals,
     * with index either a positive or negative integer to substring from left or right of the string respectively.</li>
     * If any of the indices are off-bounds, then the full chemical component identifier is replaced
     */
    public static void setChemCompPathUrlTemplate(String chemCompPathUrlTemplate) {
        DownloadChemCompProvider.chemCompPathUrlTemplate = chemCompPathUrlTemplate;
    }

    /**
     * Get this provider's cache path
     * @return
     */
    public static File getPath() {
        if (path == null) {
            UserConfiguration config = new UserConfiguration();
            path = new File(config.getCacheFilePath());
        }
        return path;
    }

    /**
     * Checks if the chemical components already have been installed into the PDB directory.
     * If not, will download the chemical components definitions file and split it up into small
     * subfiles.
     */
    public void checkDoFirstInstall() {
        if (!downloadAll) {
            return;
        }

        // this makes sure there is a file separator between every component,
        // if path has a trailing file separator or not, it will work for both cases
        File dir = new File(getPath(), CHEM_COMP_CACHE_DIRECTORY);
        File f = new File(dir, "components.cif.gz");

        if (!f.exists()) {
            downloadAllDefinitions();
        } else {
            // file exists.. did it get extracted?
            FilenameFilter filter = (dir1, file) -> file.endsWith(".cif.gz");
            String[] files = dir.list(filter);
            if (files.length < 500) {
                // not all did get unpacked
                try {
                    split();
                } catch (IOException e) {
                    logger.error("Could not split file {} into individual chemical component files. Error: {}",
                            f.toString(), e.getMessage());
                }
            }
        }
    }

    private void split() throws IOException {
        logger.info("Installing individual chem comp files ...");

        File dir = new File(getPath(), CHEM_COMP_CACHE_DIRECTORY);
        File f = new File(dir, "components.cif.gz");

        int counter = 0;
        InputStreamProvider prov = new InputStreamProvider();

        try (BufferedReader buf = new BufferedReader (new InputStreamReader(prov.getInputStream(f)))) {
            String line;
            line = buf.readLine ();
            StringWriter writer = new StringWriter();

            String currentID = null;
            while (line != null) {
                if (line.startsWith("data_")) {
                    // a new record found!

                    if (currentID != null) {
                        writeID(writer.toString(), currentID);
                        counter++;
                    }

                    currentID = line.substring(5);
                    writer = new StringWriter();
                }

                writer.append(line);
                writer.append(NEWLINE);

                line = buf.readLine();
            }

            // write the last record...
            writeID(writer.toString(), currentID);
            counter++;
        }

        logger.info("Created {} chemical component files.", counter);
    }

    /**
     * Output chemical contents to a file
     * @param contents File contents
     * @param currentID Chemical ID, used to determine the filename
     * @throws IOException
     */
    private void writeID(String contents, String currentID) throws IOException {
        String localName = getLocalFileName(currentID);
        try (PrintWriter pw = new PrintWriter(new GZIPOutputStream(new FileOutputStream(localName)))) {
            pw.print(contents);
            pw.flush();
        }
    }

    /**
     * Loads the definitions for this {@link ChemComp} from a local file and instantiates a new object.
     *
     * @param recordName the ID of the {@link ChemComp}
     * @return a new {@link ChemComp} definition.
     */
    @Override
    public ChemComp getChemComp(String recordName) {
        // make sure we work with upper case records
        recordName = recordName.toUpperCase().trim();

        boolean haveFile = true;
        if (recordName.equals("?")) {
            return null;
        }

        if (fileIsAbsent(recordName)) {
            // check if we should install all components
            checkDoFirstInstall();
        }
        if (fileIsAbsent(recordName)) {
            // we previously have installed already the definitions,
            // just do an incrememntal update
            haveFile = downloadChemCompRecord(recordName);
        }

        // Added check that download was successful and chemical component is available.
        if (haveFile) {
            String filename = getLocalFileName(recordName);
            try {
                ChemComp chemComp;
                try {
                    ChemicalComponentDictionary dict = ChemCompConverter.fromPath(Paths.get(filename));
                    chemComp = dict.getChemComp(recordName);
                } catch (ParsingException e) {
                    // happens for corrupt files
                    chemComp = null;
                }

                // May be null if the file was corrupt. Fall back on ReducedChemCompProvider in that case
                if (chemComp != null) {
                    return chemComp;
                }
            } catch (IOException e) {
                logger.warn("Could not download chemical component file {} for {}. Error: {}. Now trying to use the " +
                                "local chemical component definitions.", filename, recordName, e.getMessage());
            }
        }

        // see https://github.com/biojava/biojava/issues/315
        // probably a network error happened. Try to use the ReducedChemCOmpProvider
        if (fallback == null) {
            fallback = new ReducedChemCompProvider();
        }

        logger.warn("Falling back to ReducedChemCompProvider for {}. This could indicate a network error.", recordName);
        return fallback.getChemComp(recordName);
    }

    /**
     * Returns the file name that contains the definition for this {@link ChemComp}
     *
     * @param recordName the ID of the {@link ChemComp}
     * @return full path to the file
     */
    public static String getLocalFileName(String recordName) {
        if (protectedIDs.contains(recordName)) {
            recordName = "_" + recordName;
        }

        File f = new File(getPath(), CHEM_COMP_CACHE_DIRECTORY);
        if (!f.exists()) {
            logger.info("Creating directory {}", f);

            boolean success = f.mkdir();
            // we've checked in initPath that path is writable, so there's no need to check if it succeeds
            // in the unlikely case that in the meantime it isn't writable at least we log an error
            if (!success) {
                logger.error("Directory {} could not be created", f);
            }
        }

        File theFile = new File(f, recordName + ".cif.gz");
        return theFile.toString();
    }

    private static boolean fileIsAbsent(String recordName) {
        String fileName = getLocalFileName(recordName);
        File f = new File(fileName);

        // delete files that are too short to have contents
        if (f.length() < LocalPDBDirectory.MIN_PDB_FILE_SIZE) {
            // Delete defensively.
            // Note that if delete is unsuccessful, we re-download the file anyways
            f.delete();
            return true;
        }

        return !f.exists();
    }

    /**
     * Expands the given path URL template, replacing the placeholders as specified in {@link #setChemCompPathUrlTemplate(String)}
     * by the ccdId given (or its substrings, if indices are present in the template)
     * @param templateStr the template string with placeholders for ccd ids
     * @param ccdId the ccd id to replace (in full or a substring)
     * @return the input templateStr with placeholders replaced
     */
    static String expandPathUrlTemplate(String templateStr, String ccdId) {
        Matcher m = CCD_ID_TEMPLATE_REGEX.matcher(templateStr);
        StringBuilder output = new StringBuilder();
        int lastIndex = 0;
        while (m.find()) {
            String repString = ccdId;
            String indicesStr = m.group(1);
            try {
                if (indicesStr == null) {
                    // no substringing
                    repString = ccdId;
                } else if (!indicesStr.contains("_")) {
                    // left/right substring
                    int idx = Integer.parseInt(indicesStr);
                    if (idx < 0) { // right substring
                        repString = ccdId.substring(ccdId.length() + idx);
                    } else { // left substring
                        repString = ccdId.substring(0, idx);
                    }
                } else if (indicesStr.contains("_")) {
                    // start and end index
                    String[] tokens = indicesStr.split("_");
                    int begIdx = Integer.parseInt(tokens[0]);
                    int endIdx = Integer.parseInt(tokens[1]);
                    repString = ccdId.substring(begIdx, endIdx);
                }
            } catch (IndexOutOfBoundsException e) {
                // we don't set repString, it keeps original value ccdId
                logger.debug("Indices included in path URL template {} are out of bounds for string {}", templateStr, ccdId);
            }
            output.append(templateStr, lastIndex, m.start()).append(repString);

            lastIndex = m.end();
            // TODO when we upgrade to java 11, use the new methods introduced in java 9, see https://stackoverflow.com/questions/9605716/java-regular-expression-find-and-replace
        }
        if (lastIndex < templateStr.length()) {
            output.append(templateStr, lastIndex, templateStr.length());
        }
        return output.toString();
    }

    /**
     * @param recordName : three-letter name
     * @return true if successful download
     */
    private static boolean downloadChemCompRecord(String recordName) {
        String localName = getLocalFileName(recordName);
        File newFile;
        try {
            newFile = File.createTempFile("chemcomp" + recordName, "cif");
            logger.debug("Will write chem comp file to temp file {}", newFile.toString());
        } catch(IOException e) {
            logger.error("Could not write to temp directory {} to create the chemical component download temp file", System.getProperty("java.io.tmpdir"));
            return false;
        }

        String u = serverBaseUrl + expandPathUrlTemplate(chemCompPathUrlTemplate, recordName);

        logger.debug("Downloading chem comp definition from {}", u);

        URL url = null;
        try {
            url = new URL(u);
            URLConnection uconn = URLConnectionTools.openURLConnection(url);

            try (PrintWriter pw = new PrintWriter(new GZIPOutputStream(new FileOutputStream(newFile)));
                 BufferedReader fileBuffer = new BufferedReader(new InputStreamReader(uconn.getInputStream()))) {
                String line;
                boolean success = false;
                while ((line = fileBuffer.readLine()) != null) {
                    pw.println(line);
                    success = true;
                }
                if(!success) {
                	throw new IOException("Malformed URL or no content found in "+url.toString());
                }

                pw.flush();
            }
            // Now we move this across to where it actually wants to be
            Files.move(newFile.toPath(), Paths.get(localName), StandardCopyOption.REPLACE_EXISTING);

            return true;
        } catch (IOException e) {
            logger.error("Could not download {} OR store locally to {} Error ={}",
                    url,
                    localName,
                    e.getMessage());
            newFile.delete();
        }
        return false;
    }

    private void downloadAllDefinitions() {
        if (loading.get()) {
            logger.info("Waiting for other thread to install chemical components...");
        }

        while (loading.get()) {
            // another thread is already downloading the components definitions
            // wait for the other thread to finish...
            try {
                // wait half a second
                Thread.sleep(500);
            } catch (InterruptedException e) {
                //e.printStackTrace();
                logger.error("Thread interrupted "+e.getMessage());
            }

            logger.info("Another thread installed the chemical components.");
            return;
        }

        loading.set(true);
        long timeS = System.currentTimeMillis();

        logger.info("Performing first installation of chemical components.");
        logger.info("Downloading components.cif.gz ...");

        try {
            AllChemCompProvider.downloadFile();
        } catch (IOException e) {
            logger.error("Could not download the all chemical components file. Error: {}. "
                    + "Chemical components information won't be available", e.getMessage());
            // no point in trying to split if the file could not be downloaded
            loading.set(false);
            return;
        }
        try {
            split();
        } catch (IOException e) {
            logger.error("Could not split all chem comp file into individual chemical component files. Error: {}",
                    e.getMessage());
            // no point in reporting time
            loading.set(false);
            return;
        }
        long timeE = System.currentTimeMillis();
        logger.info("time to install chem comp dictionary: " + (timeE - timeS) / 1000 + " sec.");
        loading.set(false);
    }

    /**
     * By default this provider will download only some of the {@link ChemComp} files.
     * The user has to request that all files should be downloaded by setting this parameter to true.
     *
     *  @return flag if the all components should be downloaded and installed at startup. (default: false)
     */
    public boolean isDownloadAll() {
        return downloadAll;
    }

    /** By default this provider will download only some of the {@link ChemComp} files.
     * The user has to request that all files should be downloaded by setting this parameter to true.
     *
     * @param downloadAll if the all components should be downloaded and installed at startup. (default: false)
     */
    public void setDownloadAll(boolean downloadAll) {
        this.downloadAll = downloadAll;
    }
}
