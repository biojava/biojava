package org.biojava.nbio.structure.chem;

import org.biojava.nbio.core.util.SoftHashMap;
import org.biojava.nbio.structure.AminoAcid;
import org.biojava.nbio.structure.AminoAcidImpl;
import org.biojava.nbio.structure.Group;
import org.biojava.nbio.structure.HetatomImpl;
import org.biojava.nbio.structure.NucleotideImpl;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.Map;

public class ChemCompGroupFactory {
    private static final Logger logger = LoggerFactory.getLogger(ChemCompGroupFactory.class);
    private static ChemCompProvider chemCompProvider = new DownloadChemCompProvider();
    private static Map<String, ChemComp> cache = new SoftHashMap<>(0);

    public static ChemComp getChemComp(String recordName) {
        recordName = recordName.toUpperCase().trim();

        // we are using the cache, to avoid hitting the file system too often.
        ChemComp chemComp = cache.get(recordName);
        if (chemComp != null) {
            logger.debug("Chem comp " + chemComp.getThreeLetterCode() + " read from cache");
            return chemComp;
        }

        // not cached, get the chem comp from the provider
        logger.debug("Chem comp " + recordName + " read from provider " + chemCompProvider.getClass().getCanonicalName());
        chemComp = chemCompProvider.getChemComp(recordName);

        // Note that this also caches null or empty responses
        cache.put(recordName, chemComp);
        return chemComp;
    }

    /**
     * The new ChemCompProvider will be set in the static variable,
     * so this provider will be used from now on until it is changed
     * again. Note that this change can have unexpected behavior of
     * code executed afterwards.
     * <p>
     * Changing the provider also resets the cache, so any groups
     * previously accessed will be reread or re-downloaded.
     *
     * @param provider
     */
    public static void setChemCompProvider(ChemCompProvider provider) {
        logger.debug("Setting new chem comp provider to " + provider.getClass().getCanonicalName());
        chemCompProvider = provider;
        // clear cache
        cache.clear();
    }

    public static ChemCompProvider getChemCompProvider(){
        return chemCompProvider;
    }

    /**
     * Force the in-memory cache to be reset.
     *
     * Note that the ChemCompProvider may have additional memory or disk caches that need to be cleared too.
     */
    public static void clearCache() {
        cache.clear();
    }

    public static Group getGroupFromChemCompDictionary(String recordName) {
        // make sure we work with upper case records
        recordName = recordName.toUpperCase().trim();
        ChemComp chemComp = getChemComp(recordName);
        Group group;

        if (chemComp == null) {
            return null;
        }

        PolymerType polymerType = PolymerType.polymerTypeFromString(chemComp.getType());
        if (PolymerType.PROTEIN_ONLY.contains(polymerType)) {
            AminoAcid aminoAcid = new AminoAcidImpl();

            String oneLetterCode = chemComp.getOneLetterCode();
            if (oneLetterCode == null || oneLetterCode.equals("X") || oneLetterCode.equals("?") || oneLetterCode.length() == 0) {
                String parent = chemComp.getMonNstdParentCompId();
                if (parent != null && parent.length() == 3) {
                    String parentId = chemComp.getMonNstdParentCompId();
                    ChemComp parentChemComp = getChemComp(parentId);
                    oneLetterCode = parentChemComp.getOneLetterCode();
                }
            }

            if (oneLetterCode == null || oneLetterCode.length() == 0 || oneLetterCode.equals("?")) {
                // e.g. problem with PRR, which probably should have a parent of ALA, but as of 20110127 does not.
                logger.warn("Problem with chemical component: " + recordName + "  Did not find one letter code! Setting it to 'X'");
                aminoAcid.setAminoType('X');
            } else {
                aminoAcid.setAminoType(oneLetterCode.charAt(0));
            }

            group = aminoAcid;
        } else if (PolymerType.POLYNUCLEOTIDE_ONLY.contains(polymerType)) {
            group = new NucleotideImpl();
        } else {
            group = new HetatomImpl();
        }

        group.setChemComp(chemComp);
        return group;
    }

    public static String getOneLetterCode(ChemComp chemComp) {
        String oneLetterCode = chemComp.getOneLetterCode();
        if (oneLetterCode == null || oneLetterCode.equals("X") || oneLetterCode.equals("?")) {
            String parentId = chemComp.getMonNstdParentCompId();
            if (parentId == null) {
                return oneLetterCode;
            }
            // cases like OIM have multiple parents (comma separated), we shouldn't try grab a chemcomp for those strings
            if (parentId.length() > 3) {
                return oneLetterCode;
            }
            ChemComp parentChemComp = ChemCompGroupFactory.getChemComp(parentId);
            if (parentChemComp == null) {
                return oneLetterCode;
            }
            oneLetterCode = parentChemComp.getOneLetterCode();
        }
        return oneLetterCode;
    }
}
