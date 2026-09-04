/**
 * BioJava development code
 *
 * This code may be freely distributed and modified under the terms of the GNU
 * Lesser General Public Licence. This should be distributed with the code. If
 * you do not have a copy, see:
 *
 * http://www.gnu.org/copyleft/lesser.html
 *
 * Copyright for this code is held jointly by the individual authors. These
 * should be listed in @author doc comments.
 *
 * For more information on the BioJava project and its aims, or to join the
 * biojava-l mailing list, visit the home page at:
 *
 * http://www.biojava.org/
 */
package org.biojava.nbio.structure.io.density;

import java.io.File;

import org.biojava.nbio.structure.PdbId;
import org.biojava.nbio.structure.io.LocalPDBDirectory;

/**
 * Where cached density files live on disk.
 * <p>
 * PDB-keyed maps follow the divided layout the rest of BioJava already uses, so
 * that a cache with many entries does not end up with one enormous directory:
 * <pre>
 * &lt;cache&gt;/density/cb/1cbs_2fofc_pdbe.ccp4
 * &lt;cache&gt;/density/cb/1cbs_fofc_pdbe.ccp4
 * &lt;cache&gt;/density/cb/1cbs_2fofc_rcsbvs_d0.bcif
 * &lt;cache&gt;/density/cb/1cbs_2fofc_wwpdb.cif.gz
 * </pre>
 * Both the kind and the source appear in the file name, so no two combinations
 * can collide.
 * <p>
 * EMDB maps are keyed by EMDB identifier instead, mirroring the EMDB archive:
 * <pre>
 * &lt;cache&gt;/density/emd/EMD-0262/emd_0262.map.gz
 * </pre>
 * Several PDB entries are often fitted into a single EM map, and those maps can
 * be hundreds of megabytes, so keying them by PDB entry would cache the same
 * enormous file many times over.
 * <p>
 * The two-character directory needs no attention when the archive moves to
 * extended identifiers in July 2027. It is a device for spreading files over
 * directories, not a copy of the archive's own layout, and
 * {@link LocalPDBDirectory#getMiddleHash(String)} counts from the right hand end
 * of the identifier, so <code>1cbs</code> and <code>pdb_00001cbs</code> both land
 * in <code>cb</code> &mdash; which is also the rule the wwPDB documents for the
 * new archive.
 * <p>
 * Deliberately <i>not</i> the archive's per-entry layout. A cache is not a mirror
 * and cannot become one: the archive publishes structure factors and map
 * coefficients but never grids, so density has to be fetched whatever else is
 * mirrored, and writing it into a directory whose contents are an exact copy of
 * upstream puts it at the mercy of the next <code>rsync --delete</code>. Should
 * that judgement ever be revisited, every cached path is computed here and
 * nowhere else, so a different layout is a change to this class plus a fallback
 * probe for files in the old places &mdash; not a cache that everyone has to
 * discard.
 *
 * @author Amr ALHOSSARY
 * @since 7.3.0
 */
public class DensityCacheLayout {

	/** Name of the density sub-directory within the BioJava cache directory. */
	public static final String DENSITY_DIR = "density";

	/** Name of the sub-directory holding EMDB-keyed maps. */
	public static final String EMDB_DIR = "emd";

	/** Name of the sub-directory holding cached PDB-to-EMDB mappings. */
	public static final String EMDB_MAPPING_DIR = "emdb-mapping";

	/**
	 * File-name token for a file that holds more than one kind of map, as a density
	 * server response does.
	 */
	public static final String BOTH_KINDS_TOKEN = "both";

	private DensityCacheLayout() {
	}

	/**
	 * The root density directory inside a cache directory.
	 *
	 * @param cacheRoot the BioJava cache directory
	 * @return the density directory, which need not exist yet
	 */
	public static File densityRoot(File cacheRoot) {
		return new File(cacheRoot, DENSITY_DIR);
	}

	/**
	 * The cache file for a PDB-keyed map.
	 *
	 * @param cacheRoot the BioJava cache directory
	 * @param pdbId the entry
	 * @param kind the kind of map
	 * @param source the service it came from
	 * @param format the file format
	 * @param qualifier an extra discriminator such as a detail level, or
	 *        <code>null</code>. Anything that changes the content but not the
	 *        entry, kind or source belongs here.
	 * @return the file, which need not exist
	 */
	public static File pdbMapFile(File cacheRoot, PdbId pdbId, DensityMapKind kind, DensityMapSource source,
			DensityFileFormat format, String qualifier) {
		return pdbMapFile(cacheRoot, pdbId, kind.getFileToken(), source, format, qualifier);
	}

	/**
	 * The cache file for a PDB-keyed map, naming the kind explicitly.
	 * <p>
	 * The separate kind token exists for sources that deliver more than one kind of
	 * map in a single file. A density server response, for instance, carries both
	 * the 2Fo-Fc and the Fo-Fc blocks, so it is cached once under
	 * {@link #BOTH_KINDS_TOKEN} rather than downloaded and stored twice.
	 *
	 * @param cacheRoot the BioJava cache directory
	 * @param pdbId the entry
	 * @param kindToken the token naming what the file holds
	 * @param source the service it came from
	 * @param format the file format
	 * @param qualifier an extra discriminator such as a detail level, or <code>null</code>
	 * @return the file, which need not exist
	 */
	public static File pdbMapFile(File cacheRoot, PdbId pdbId, String kindToken, DensityMapSource source,
			DensityFileFormat format, String qualifier) {
		String id = shortIdOrFull(pdbId).toLowerCase();
		File dir = new File(densityRoot(cacheRoot), LocalPDBDirectory.getMiddleHash(id));
		StringBuilder name = new StringBuilder(id)
				.append('_').append(kindToken)
				.append('_').append(source.getFileToken());
		if (qualifier != null && !qualifier.isEmpty()) {
			name.append('_').append(qualifier);
		}
		name.append(format.getExtension());
		return new File(dir, name.toString());
	}

	/**
	 * The cache file for an EMDB-keyed map.
	 *
	 * @param cacheRoot the BioJava cache directory
	 * @param emdbId the EMDB entry, in any accepted form
	 * @param source the service it came from
	 * @param format the file format
	 * @param qualifier an extra discriminator such as a detail level, or <code>null</code>
	 * @return the file, which need not exist
	 */
	public static File emdbMapFile(File cacheRoot, String emdbId, DensityMapSource source,
			DensityFileFormat format, String qualifier) {
		String canonical = DensityMapRequest.normalizeEmdbId(emdbId);
		String number = DensityMapRequest.emdbNumber(emdbId);
		File dir = new File(new File(densityRoot(cacheRoot), EMDB_DIR), canonical);
		StringBuilder name = new StringBuilder("emd_").append(number);
		if (source != DensityMapSource.EMDB_MAP) {
			name.append('_').append(source.getFileToken());
		}
		if (qualifier != null && !qualifier.isEmpty()) {
			name.append('_').append(qualifier);
		}
		name.append(format.getExtension());
		return new File(dir, name.toString());
	}

	/**
	 * The file caching the EMDB identifiers and author contour level associated
	 * with a PDB entry.
	 *
	 * @param cacheRoot the BioJava cache directory
	 * @param pdbId the entry
	 * @return the file, which need not exist
	 */
	public static File emdbMappingFile(File cacheRoot, PdbId pdbId) {
		String id = shortIdOrFull(pdbId).toLowerCase();
		File dir = new File(new File(densityRoot(cacheRoot), EMDB_MAPPING_DIR), LocalPDBDirectory.getMiddleHash(id));
		return new File(dir, id + ".emdb.properties");
	}

	/**
	 * The file caching an EMDB entry's map metadata as served by the EMDB API.
	 *
	 * @param cacheRoot the BioJava cache directory
	 * @param emdbId the EMDB entry, in any accepted form
	 * @return the file, which need not exist
	 */
	public static File emdbMapInfoFile(File cacheRoot, String emdbId) {
		String canonical = DensityMapRequest.normalizeEmdbId(emdbId);
		File dir = new File(new File(densityRoot(cacheRoot), EMDB_DIR), canonical);
		return new File(dir, canonical + ".map-info.json");
	}

	/**
	 * The companion name a density-server file must have for Jmol to read its
	 * difference-map block.
	 * <p>
	 * A density server response holds both a <code>2FO-FC</code> and an
	 * <code>FO-FC</code> data block, and Jmol's reader chooses between them by
	 * testing whether the <i>file name</i> contains the literal text
	 * <code>&amp;diff=1</code>. That works for a URL fetched straight from the
	 * server, where the marker rides along in the query string, but a cached local
	 * file has no query string; appending the marker to the file URL only makes
	 * Jmol look for a file that does not exist. Putting the marker into the name
	 * itself is what actually selects the block.
	 * <p>
	 * This was verified against Jmol 14.31.10 and is unchanged in current Jmol: the
	 * relevant line in <code>BCifDensityReader</code> still carries the author's
	 * "what about cached data" to-do beside it. Should Jmol gain a cleaner way to
	 * choose the block, this can be retired.
	 *
	 * @param mapFile the cached density-server file
	 * @return the sibling path that selects the difference map
	 */
	public static File differenceMarkerFile(File mapFile) {
		String name = mapFile.getName();
		int dot = name.lastIndexOf('.');
		String stem = dot < 0 ? name : name.substring(0, dot);
		String ext = dot < 0 ? "" : name.substring(dot);
		return new File(mapFile.getAbsoluteFile().getParentFile(), stem + "&diff=1" + ext);
	}

	/**
	 * The short four-character spelling of an identifier where one exists.
	 * <p>
	 * Every entry in the archive today has a four-character form, and the density
	 * services accept only that spelling. Extended-only identifiers will appear
	 * eventually; rather than refusing them, the extended spelling is passed
	 * through so that the services can start accepting it without a change here.
	 *
	 * @param pdbId the identifier
	 * @return the short spelling if available, otherwise the full one
	 */
	public static String shortIdOrFull(PdbId pdbId) {
		return pdbId.getId(true);
	}
}
