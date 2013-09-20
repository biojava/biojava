package org.biojava.bio.structure.scop;




/**
 * Controls the global ScopDatabase being used.
 *
 * <p>Defaults to a {@link RemoteScopInstallation}, which is fast for small numbers
 * of queries. For many queries, a {@link ScopInstallation} instance may be faster,
 * since it makes only one network request.
 *
 * <p>Example: Fetch the structure corresponding to an old version of scop
 *
 * <pre>
 * ScopInstallation scop = new ScopInstallation();
 * scop.setScopVersion("1.69");
 * ScopFactory.setScopDatabase(scop);
 * AtomCache cache = new AtomCache();
 * cache.setFetchFileEvenIfObsolete(true); //fetch older PDBs
 * cache.setStrictSCOP(false); // correct simple errors in domain names
 * Structure s = cache.getStructure("d3hbia_");
 * @author sbliven
 *
 */
public class ScopFactory {

    // berkeley
    public static final String VERSION_1_75B = "1.75B";
    public static final String VERSION_1_75A = "1.75A";
    // original SCOP
    // latest SCOP release from SCOP website = 1.75;
    public static final String VERSION_1_75 = "1.75";
    // outdated stable versions
    public static final String VERSION_1_73 = "1.73";
    public static final String VERSION_1_71 = "1.71";
    public static final String VERSION_1_69 = "1.69";
    public static final String VERSION_1_67 = "1.67";
    public static final String VERSION_1_65 = "1.65";
    public static final String VERSION_1_63 = "1.63";
    public static final String VERSION_1_61 = "1.61";
    public static final String VERSION_1_59 = "1.59";
    public static final String VERSION_1_57 = "1.57";
    public static final String VERSION_1_55 = "1.55";

    // by default we now request SCOP via web services -> much less memory consumption
    private static ScopDatabase scop = new RemoteScopInstallation();

    /**
     * Get the current default instance
     * @return
     */
    public static ScopDatabase getSCOP(){
        return scop;
    }

    /**
     * requests a particular version of SCOP.
     *
     * Where possible, this will be the current default instance.
     * Otherwise a new instance will be created.
     * @param version
     * @return
     */
    public static ScopDatabase getSCOP(String version){
        // assume versions sort lexicographically
        if(scop.getScopVersion().equalsIgnoreCase(version)) {
            return scop;
        }
        // Default to a local installation
        return getSCOP(version,true);
    }

    /**
     * Gets an instance of the specified scop version. Not that the instance
     * may or may not be the default instance.
     *
     * If useLocalData, the returned instance will be either a
     * {@link ScopInstallation} or a {@link BerkeleyScopInstallation}. If not,
     * it will be a {@link RemoteScopInstallation}.
     * @param version A version number, such as {@link #VERSION_1_75A}
     * @param useLocalData Whether to use a local installation or a remote installation
     * @return an
     */
    public static ScopDatabase getSCOP(String version, boolean useLocalData){
        if ( useLocalData) {
            // Use a local installation
            // Assume version strings sort lexicographically
            if ( version.compareToIgnoreCase(VERSION_1_75) > 0 ) {
                return getBerkeley(version);
            } else {
                return getScopInstallation(version);
            }
        } else {
            // Use a remote installation
            return scop;
        }
    }

    /**
     * Set the default instance to use the specified SCOP version number
     * @param version A version number, such as {@link #VERSION_1_75A}
     * @return the new default instance
     */
    public static ScopDatabase setScopDatabase(String version) {
        setScopDatabase(ScopFactory.getSCOP(version));
        return scop;
    }

    /**
     * Make `scop` the default instance for all BioJava
     * @param scop
     */
    public static void setScopDatabase(ScopDatabase scop){
        System.out.println("Setting ScopDatabase to type: " + scop.getClass().getName());
        ScopFactory.scop = scop;
    }

    /**
     * Gets a local instance with the specified version.
     *
     * Uses the singleton if applicable, or creates a new instance
     * @param version a version applicable to a ScopInstallation
     * @return
     */
    private static ScopInstallation getScopInstallation(String version) {
        if( scop instanceof ScopInstallation && scop.getScopVersion() == version) {
            return (ScopInstallation)scop;
        } else {
            ScopInstallation scopInst = new ScopInstallation();
            scopInst.setScopVersion(version);
            return scopInst;
        }
    }

    /**
     * Gets a Berkeley instance with the specified version.
     *
     * Uses the singleton if applicable, or creates a new instance
     * @param version
     * @return
     */
    private static BerkeleyScopInstallation getBerkeley(String version) {
        if( scop instanceof BerkeleyScopInstallation && scop.getScopVersion() == version) {
            return (BerkeleyScopInstallation)scop;
        } else {
            BerkeleyScopInstallation berkeley  = new BerkeleyScopInstallation();
            berkeley.setScopVersion(version);
            return berkeley;
        }
    }
}
