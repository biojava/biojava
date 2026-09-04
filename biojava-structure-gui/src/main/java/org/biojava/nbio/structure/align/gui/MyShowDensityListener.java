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
package org.biojava.nbio.structure.align.gui;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ExecutionException;

import javax.swing.JOptionPane;
import javax.swing.SwingWorker;

import org.biojava.nbio.structure.PdbId;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.align.gui.jmol.AbstractAlignmentJmol;
import org.biojava.nbio.structure.io.density.DensityMapCache;
import org.biojava.nbio.structure.io.density.DensityMapKind;
import org.biojava.nbio.structure.io.density.DensityMapRequest;
import org.biojava.nbio.structure.io.density.DensityMapResult;
import org.biojava.nbio.structure.io.density.DensityMapSource;
import org.biojava.nbio.structure.io.density.NoDensityMapException;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Fetches the density map for the structure on screen and draws it.
 * <p>
 * The download happens on a {@link SwingWorker} rather than the event dispatch
 * thread. Even the smallest source runs to a few hundred kilobytes and a
 * full-resolution map can be far larger, so fetching inline would freeze the
 * interface for as long as it took.
 *
 * @author Amr ALHOSSARY
 * @since 7.3.0
 */
public class MyShowDensityListener implements ActionListener {

	private static final Logger logger = LoggerFactory.getLogger(MyShowDensityListener.class);

	private static final String OPTION_2FOFC = "2Fo-Fc (electron density)";
	private static final String OPTION_FOFC = "Fo-Fc (difference)";
	private static final String OPTION_BOTH = "Both";
	private static final String OPTION_AUTO = "Whatever is available";

	private final AbstractAlignmentJmol parent;

	/**
	 * @param parent the viewer to draw into
	 */
	public MyShowDensityListener(AbstractAlignmentJmol parent) {
		this.parent = parent;
	}

	@Override
	public void actionPerformed(ActionEvent e) {
		Structure structure = parent == null ? null : parent.getStructure();
		if (structure == null) {
			JOptionPane.showMessageDialog(parent == null ? null : parent.getFrame(),
					"There is no structure on screen to fetch a density map for.",
					"Show Electron Density", JOptionPane.INFORMATION_MESSAGE);
			return;
		}

		PdbId pdbId = structure.getPdbId();
		if (pdbId == null) {
			String typed = JOptionPane.showInputDialog(parent.getFrame(),
					"This structure has no PDB ID. Enter one to look up its density map:",
					"Show Electron Density", JOptionPane.QUESTION_MESSAGE);
			if (typed == null || typed.trim().isEmpty()) {
				return;
			}
			try {
				pdbId = new PdbId(typed.trim());
			} catch (IllegalArgumentException ex) {
				JOptionPane.showMessageDialog(parent.getFrame(), typed + " is not a valid PDB ID.",
						"Show Electron Density", JOptionPane.ERROR_MESSAGE);
				return;
			}
		}

		Object[] options = {OPTION_2FOFC, OPTION_FOFC, OPTION_BOTH, OPTION_AUTO};
		Object choice = JOptionPane.showInputDialog(parent.getFrame(),
				"Which map would you like to see for " + pdbId.getId() + "?",
				"Show Electron Density", JOptionPane.QUESTION_MESSAGE, null, options, OPTION_2FOFC);
		if (choice == null) {
			return;
		}

		List<DensityMapKind> kinds = new ArrayList<>(2);
		if (OPTION_BOTH.equals(choice)) {
			kinds.add(DensityMapKind.TWO_FO_FC);
			kinds.add(DensityMapKind.FO_FC);
		} else if (OPTION_FOFC.equals(choice)) {
			kinds.add(DensityMapKind.FO_FC);
		} else if (OPTION_AUTO.equals(choice)) {
			kinds.add(DensityMapKind.AUTO);
		} else {
			kinds.add(DensityMapKind.TWO_FO_FC);
		}

		fetchAndShow(pdbId, kinds);
	}

	private void fetchAndShow(PdbId pdbId, List<DensityMapKind> kinds) {
		parent.setStatus("Fetching density map for " + pdbId.getId() + " ...");

		new SwingWorker<List<DensityMapResult>, Void>() {

			private NoDensityMapException missing;

			@Override
			protected List<DensityMapResult> doInBackground() throws Exception {
				DensityMapCache cache = DensityMapCache.getInstance();
				List<DensityMapResult> results = new ArrayList<>(kinds.size());
				for (DensityMapKind kind : kinds) {
					try {
						// Restricting to displayable formats also excludes the map
						// coefficient source automatically.
						results.add(cache.getDensityMap(DensityMapRequest.builder(pdbId)
								.kind(kind)
								.allowNonRenderableFormats(false)
								.build()));
					} catch (NoDensityMapException ex) {
						missing = ex;
					}
				}
				return results;
			}

			@Override
			protected void done() {
				List<DensityMapResult> results;
				try {
					results = get();
				} catch (InterruptedException ex) {
					Thread.currentThread().interrupt();
					return;
				} catch (ExecutionException ex) {
					logger.error("Could not fetch a density map for {}", pdbId.getId(), ex.getCause());
					parent.setStatus("Could not fetch density map");
					JOptionPane.showMessageDialog(parent.getFrame(),
							"Could not fetch the density map:\n" + ex.getCause().getMessage(),
							"Show Electron Density", JOptionPane.ERROR_MESSAGE);
					return;
				}

				for (DensityMapResult result : results) {
					parent.getJmolPanel().loadDensityMap(result);
				}

				if (results.isEmpty()) {
					parent.setStatus("No density map available");
					JOptionPane.showMessageDialog(parent.getFrame(), explain(pdbId, missing),
							"Show Electron Density", JOptionPane.INFORMATION_MESSAGE);
				} else {
					DensityMapResult first = results.get(0);
					parent.setStatus(String.format("Density from %s (%,d kB)",
							first.getSource(), first.getFileSizeBytes() / 1024));
				}
			}
		}.execute();
	}

	/**
	 * Turns the per-source reasons into something a user can act on, rather than a
	 * list of HTTP codes.
	 */
	private static String explain(PdbId pdbId, NoDensityMapException missing) {
		StringBuilder sb = new StringBuilder("No density map is available for ")
				.append(pdbId.getId()).append(".\n\n");
		if (missing == null) {
			return sb.toString();
		}
		Map<DensityMapSource, String> attempts = missing.getAttempts();
		boolean allNotFound = !attempts.isEmpty() && attempts.values().stream()
				.allMatch(reason -> reason.startsWith("HTTP 404") || reason.startsWith("no "));
		if (allNotFound) {
			sb.append("The most likely reason is that no structure factors were deposited\n")
			  .append("for this entry, and that it has no associated EMDB map.\n\n");
		}
		sb.append("Sources tried:\n");
		attempts.forEach((source, reason) -> sb.append("    ").append(source).append(": ").append(reason).append('\n'));
		return sb.toString();
	}
}
