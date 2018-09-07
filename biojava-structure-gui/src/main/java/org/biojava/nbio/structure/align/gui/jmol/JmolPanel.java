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
 * Created on Oct 6, 2009
 * Author: Andreas Prlic
 *
 */

package org.biojava.nbio.structure.align.gui.jmol;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.Rectangle;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.BufferedInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.PipedInputStream;
import java.io.PipedOutputStream;
import java.text.DecimalFormat;
import java.util.List;

import javax.swing.JComboBox;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Calc;
import org.biojava.nbio.structure.Group;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.align.gui.JPrintPanel;
import org.biojava.nbio.structure.domain.LocalProteinDomainParser;
import org.biojava.nbio.structure.domain.pdp.Domain;
import org.biojava.nbio.structure.domain.pdp.Segment;
import org.biojava.nbio.structure.gui.util.color.ColorUtils;
import org.biojava.nbio.structure.io.mmtf.MmtfActions;
import org.biojava.nbio.structure.jama.Matrix;
import org.biojava.nbio.structure.scop.ScopDatabase;
import org.biojava.nbio.structure.scop.ScopDomain;
import org.biojava.nbio.structure.scop.ScopInstallation;
import org.jmol.adapter.smarter.SmarterJmolAdapter;
import org.jmol.api.JmolAdapter;
import org.jmol.api.JmolStatusListener;
import org.jmol.api.JmolViewer;
import org.jmol.util.LoggerInterface;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;


public class JmolPanel
extends JPrintPanel
implements ActionListener
{
	private static final Logger logger = LoggerFactory.getLogger(JmolPanel.class);
	
	private static final long serialVersionUID = -3661941083797644242L;

	private JmolViewer viewer;
	private JmolAdapter adapter;
	JmolStatusListener statusListener;
	final Dimension currentSize = new Dimension();
	final Rectangle rectClip = new Rectangle();

	Structure structure;

	public JmolPanel() {
		super();
		statusListener = new MyJmolStatusListener();
		adapter = new SmarterJmolAdapter();
		JmolLoggerAdapter jmolLogger = new JmolLoggerAdapter(LoggerFactory.getLogger(org.jmol.api.JmolViewer.class));
		org.jmol.util.Logger.setLogger(jmolLogger);
		org.jmol.util.Logger.setLogLevel( jmolLogger.getLogLevel() );
		viewer = JmolViewer.allocateViewer(this,
				adapter,
				null,null,null,null,
				statusListener);

	}

	@Override
	public void paint(Graphics g) {

		getSize(currentSize);
		g.getClipBounds(rectClip);
		viewer.renderScreenImage(g, currentSize, rectClip);
	}

	public void evalString(String rasmolScript){

		viewer.evalString(rasmolScript);

	}

	public void openStringInline(String pdbFile){
		viewer.openStringInline(pdbFile);

	}
	public JmolViewer getViewer() {
		return viewer;
	}

	public JmolAdapter getAdapter(){
		return adapter;
	}

	public JmolStatusListener getStatusListener(){
		return statusListener;
	}
	public void executeCmd(String rasmolScript) {
		viewer.evalString(rasmolScript);
	}
	
	public void setStructure(final Structure s, boolean useMmtf) {
		
		this.structure = s;
		
		if (useMmtf) {
			try (
					PipedOutputStream out = new PipedOutputStream();
					// Viewer requires a BufferedInputStream for reflection
					InputStream in = new BufferedInputStream(new PipedInputStream(out));
					) {
				new Thread((Runnable)() -> {
					try {
						MmtfActions.writeToOutputStream(s,out);
					} catch (Exception e) {
						logger.error("Error generating MMTF output for {}",
								s.getStructureIdentifier()==null ? s.getStructureIdentifier().getIdentifier() : s.getName(), e);
					}
				}).start();
				viewer.openReader(null, in);
			} catch (IOException e) {
				logger.error("Error transfering {} to Jmol",
						s.getStructureIdentifier()==null ? s.getStructureIdentifier().getIdentifier() : s.getName(), e);
			}
		} else {
			// Use mmCIF format
			String serialized = s.toMMCIF();
			viewer.openStringInline(serialized);
			
		}
		
		evalString("save STATE state_1");
		
	}

	public void setStructure(final Structure s) {
		// Set the default to MMCIF (until issue #629 is fixed)
		setStructure(s, false);
	}

	/** assign a custom color to the Jmol chains command.
	 *
	 */
	public void jmolColorByChain(){
		String script =
				"function color_by_chain(objtype, color_list) {"+ String.format("%n") +
				""+ String.format("%n") +
				"		 if (color_list) {"+ String.format("%n") +
				"		   if (color_list.type == \"string\") {"+ String.format("%n") +
				"		     color_list = color_list.split(\",\").trim();"+ String.format("%n") +
				"		   }"+ String.format("%n") +
				"		 } else {"+ String.format("%n") +
				"		   color_list = [\"104BA9\",\"AA00A2\",\"C9F600\",\"FFA200\",\"284A7E\",\"7F207B\",\"9FB82E\",\"BF8B30\",\"052D6E\",\"6E0069\",\"83A000\",\"A66A00\",\"447BD4\",\"D435CD\",\"D8FA3F\",\"FFBA40\",\"6A93D4\",\"D460CF\",\"E1FA71\",\"FFCC73\"];"+ String.format("%n") +
				"		 }"+ String.format("%n") +

				"		 var cmd2 = \"\";"+ String.format("%n") +

				"		 if (!objtype) {"+ String.format("%n") +
				"		   var type_list  = [ \"backbone\",\"cartoon\",\"dots\",\"halo\",\"label\",\"meshribbon\",\"polyhedra\",\"rocket\",\"star\",\"strand\",\"strut\",\"trace\"];"+ String.format("%n") +
				"		   cmd2 = \"color \" + type_list.join(\" none; color \") + \" none;\";"+ String.format("%n") +
				"		   objtype = \"atoms\";"+ String.format("%n") +

				"		 }"+ String.format("%n") +

				"		 var chain_list  = script(\"show chain\").trim().lines;"+ String.format("%n") +
				"		 var chain_count = chain_list.length;"+ String.format("%n") +

				"		 var color_count = color_list.length;"+ String.format("%n") +
				"		 var sel = {selected};"+ String.format("%n") +
				"		 var cmds = \"\";"+ String.format("%n") +


				"		 for (var chain_number=1; chain_number<=chain_count; chain_number++) {"+ String.format("%n") +
				"		   // remember, Jmol arrays start with 1, but % can return 0"+ String.format("%n") +
				"		   cmds += \"select sel and :\" + chain_list[chain_number] + \";color \" + objtype + \" [x\" + color_list[(chain_number-1) % color_count + 1] + \"];\" + cmd2;"+ String.format("%n") +
				"		 }"+ String.format("%n") +
				"		 script INLINE @{cmds + \"select sel\"}"+ String.format("%n") +
				"}";

		executeCmd(script);
	}

	/** The user selected one of the Combo boxes...
	 *
	 * @param event an ActionEvent
	 */
	@Override
	public void actionPerformed(ActionEvent event) {

		Object mysource = event.getSource();

		if ( ! (mysource instanceof JComboBox )) {
			super.actionPerformed(event);
			return;
		}

		@SuppressWarnings("unchecked")
		JComboBox<String> source = (JComboBox<String>) event.getSource();
		String value = source.getSelectedItem().toString();
		evalString("save selection; ");

		String selectLigand = "select ligand;wireframe 0.16;spacefill 0.5; color cpk ;";

		if ( value.equals("Cartoon")){
			String script = "hide null; select all;  spacefill off; wireframe off; backbone off;" +
					" cartoon on; " +
					" select ligand; wireframe 0.16;spacefill 0.5; color cpk; " +
					" select *.FE; spacefill 0.7; color cpk ; " +
					" select *.CU; spacefill 0.7; color cpk ; " +
					" select *.ZN; spacefill 0.7; color cpk ; " +
					" select all; ";
			this.executeCmd(script);
		} else if (value.equals("Backbone")){
			String script = "hide null; select all; spacefill off; wireframe off; backbone 0.4;" +
					" cartoon off; " +
					" select ligand; wireframe 0.16;spacefill 0.5; color cpk; " +
					" select *.FE; spacefill 0.7; color cpk ; " +
					" select *.CU; spacefill 0.7; color cpk ; " +
					" select *.ZN; spacefill 0.7; color cpk ; " +
					" select all; ";
			this.executeCmd(script);
		} else if (value.equals("CPK")){
			String script = "hide null; select all; spacefill off; wireframe off; backbone off;" +
					" cartoon off; cpk on;" +
					" select ligand; wireframe 0.16;spacefill 0.5; color cpk; " +
					" select *.FE; spacefill 0.7; color cpk ; " +
					" select *.CU; spacefill 0.7; color cpk ; " +
					" select *.ZN; spacefill 0.7; color cpk ; " +
					" select all; ";
			this.executeCmd(script);

		} else if (value.equals("Ligands")){
			this.executeCmd("restrict ligand; cartoon off; wireframe on;  display selected;");
		} else if (value.equals("Ligands and Pocket")){
			this.executeCmd(" select within (6.0,ligand); cartoon off; wireframe on; backbone off; display selected; ");
		} else if ( value.equals("Ball and Stick")){
			String script = "hide null; restrict not water;  wireframe 0.2; spacefill 25%;" +
					" cartoon off; backbone off; " +
					" select ligand; wireframe 0.16; spacefill 0.5; color cpk; " +
					" select *.FE; spacefill 0.7; color cpk ; " +
					" select *.CU; spacefill 0.7; color cpk ; " +
					" select *.ZN; spacefill 0.7; color cpk ; " +
					" select all; ";
			this.executeCmd(script);
		} else if ( value.equals("By Chain")){
			jmolColorByChain();
			String script = "hide null; select all;set defaultColors Jmol; color_by_chain(\"cartoon\"); color_by_chain(\"\"); " + selectLigand + "; select all; ";
			this.executeCmd(script);
		} else if ( value.equals("Rainbow")) {
			this.executeCmd("hide null; select all; set defaultColors Jmol; color group; color cartoon group; " + selectLigand + "; select all; " );
		} else if ( value.equals("Secondary Structure")){
			this.executeCmd("hide null; select all; set defaultColors Jmol; color structure; color cartoon structure;" + selectLigand + "; select all; " );

		} else if ( value.equals("By Element")){
			this.executeCmd("hide null; select all; set defaultColors Jmol; color cpk; color cartoon cpk; " + selectLigand + "; select all; ");
		} else if ( value.equals("By Amino Acid")){
			this.executeCmd("hide null; select all; set defaultColors Jmol; color amino; color cartoon amino; " + selectLigand + "; select all; " );
		} else if ( value.equals("Hydrophobicity") ){
			this.executeCmd("hide null; set defaultColors Jmol; select hydrophobic; color red; color cartoon red; select not hydrophobic ; color blue ; color cartoon blue; "+ selectLigand+"; select all; ");
		} else if ( value.equals("Suggest Domains")){
			colorByPDP();
		} else if ( value.equals("Show SCOP Domains")){
			colorBySCOP();
		}
		evalString("restore selection; ");
	}

	private void colorBySCOP() {

		if ( structure == null)
			return;

		String pdbId = structure.getPDBCode();
		if ( pdbId == null)
			return;
		ScopDatabase scop = new ScopInstallation();

		List<ScopDomain> domains = scop.getDomainsForPDB(pdbId);
		if ( domains == null) {
			System.err.println("No SCOP domains found for " + pdbId);
			return;
		}
		int i = -1;
		for ( ScopDomain domain : domains){
			i++;
			if ( i >= ColorUtils.colorWheel.length)
				i = 0;
			Color c1 = ColorUtils.colorWheel[i];
			List<String>ranges = domain.getRanges();

			for (String range : ranges){
				logger.debug(range);
				String[] spl = range.split(":");
				String script = " select  ";
				if ( spl.length > 1 )
					script += spl[1]+":"+spl[0] +"/1;";
				else
					script += "*" + spl[0]+"/1;";
				script += " color [" + c1.getRed() + ","+c1.getGreen() + "," +c1.getBlue()+"];";
				script += " color cartoon [" + c1.getRed() + ","+c1.getGreen() + "," +c1.getBlue()+"] ;";
				logger.debug(script);
				evalString(script);

			}
		}


	}

	private void colorByPDP() {
		logger.debug("colorByPDP");
		if ( structure == null)
			return;

		try {
			Atom[] ca = StructureTools.getRepresentativeAtomArray(structure);
			List<Domain> domains = LocalProteinDomainParser.suggestDomains(ca);
			int i = -1;
			for ( Domain dom : domains){
				i++;
				if ( i > ColorUtils.colorWheel.length)
					i = 0;
				//System.out.println("DOMAIN:" + i + " size:" + dom.size + " " +  dom.score);
				List<Segment> segments = dom.getSegments();
				Color c1 = ColorUtils.colorWheel[i];
				//float fraction = 0f;
				for ( Segment s : segments){
					//System.out.println("   Segment: " + s);
					//Color c1 = ColorUtils.rotateHue(c, fraction);
					//	fraction += 1.0/(float)segments.size();
					int start = s.getFrom();
					int end = s.getTo();
					Group startG = ca[start].getGroup();
					Group endG = ca[end].getGroup();
					logger.debug("   Segment: " +startG.getResidueNumber() +":" + startG.getChainId() + " - " + endG.getResidueNumber()+":"+endG.getChainId() + " " + s);
					String j1 = startG.getResidueNumber()+"";
					String j2 = endG.getResidueNumber()+":"+endG.getChainId();
					String script = " select  " +j1 +"-" +j2 +"/1;";
					script += " color [" + c1.getRed() + ","+c1.getGreen() + "," +c1.getBlue()+"];";
					script += " color cartoon [" + c1.getRed() + ","+c1.getGreen() + "," +c1.getBlue()+"] ;";
					logger.debug(script);
					evalString(script);
				}

			}
		} catch (Exception e){
			e.printStackTrace();
		}
	}

	public void rotateJmol(Matrix jmolRotation) {

		if ( jmolRotation != null) {
			double[] zyz = Calc.getZYZEuler(jmolRotation);
			DecimalFormat df = new DecimalFormat("0.##");

			String script = "reset; rotate z "
					+ df.format(zyz[0])
					+ "; rotate y "
					+ df.format(zyz[1])
					+"; rotate z "
					+ df.format(zyz[2])+";";

			executeCmd(script);

		}
	}

	/** Clean up this instance for garbage collection, to avoid memory leaks...
	 *
	 */
	public void destroy(){

		executeCmd("zap;");
		structure = null;

		viewer = null;
		adapter = null;
	}

	public static class JmolLoggerAdapter implements LoggerInterface {
		private Logger slf;
		public JmolLoggerAdapter(Logger slf) {
			this.slf=slf;
		}
		public int getLogLevel() {
			if( slf.isTraceEnabled() )
				return org.jmol.util.Logger.LEVEL_MAX;
			if( slf.isDebugEnabled() )
				return org.jmol.util.Logger.LEVEL_DEBUG;
			if( slf.isInfoEnabled() )
				return org.jmol.util.Logger.LEVEL_INFO;
			if( slf.isWarnEnabled() )
				return org.jmol.util.Logger.LEVEL_WARN;
			if( slf.isErrorEnabled() )
				return org.jmol.util.Logger.LEVEL_ERROR;
			throw new IllegalStateException("Unknown SLF4J error level");
		}
		@Override
		public void debug(String txt) {
			slf.debug(txt);
		}
		@Override
		public void info(String txt) {
			slf.info(txt);
		}
		@Override
		public void warn(String txt) {
			slf.warn(txt);
		}
		@Override
		public void warnEx(String txt, Throwable e) {
			slf.warn(txt,e);
		}
		@Override
		public void error(String txt) {
			slf.error(txt);
		}
		@Override
		public void errorEx(String txt, Throwable e) {
			slf.error(txt,e);
		}
		@Override
		public void fatal(String txt) {
			slf.error(txt);
		}
		@Override
		public void fatalEx(String txt, Throwable e) {
			slf.error(txt,e);
		}
	}
}
