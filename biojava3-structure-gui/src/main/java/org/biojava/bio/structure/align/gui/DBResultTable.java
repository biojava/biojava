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
 * Created on Nov 6, 2009
 * Author: Andreas Prlic 
 *
 */

package org.biojava.bio.structure.align.gui;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.URL;
import java.util.ArrayList;
import java.util.List;


import javax.swing.JFrame;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.ListSelectionModel;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
import javax.swing.table.TableModel;
import javax.swing.table.TableRowSorter;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.StructureTools;
import org.biojava.bio.structure.align.StructureAlignment;
import org.biojava.bio.structure.align.StructureAlignmentFactory;
import org.biojava.bio.structure.align.ce.CeMain;
import org.biojava.bio.structure.align.ce.CeParameters;
import org.biojava.bio.structure.align.ce.ConfigStrucAligParams;
import org.biojava.bio.structure.align.gui.jmol.StructureAlignmentJmol;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.align.util.UserConfiguration;
import org.biojava.bio.structure.align.webstart.WebStartMain;
import org.biojava.bio.structure.io.PDBFileReader;
import org.biojava.bio.structure.io.StructureIOFile;


public class DBResultTable implements ActionListener{

	public static final String[] ceColumnNames =  {"name1","tname2","score","z-score"    ,"rmsd","len1","len2","cov1","cov2","%ID","Description",""};
	public static final String[] fatColumnNames = {"name1","tname2","score","probability","rmsd","len1","len2","cov1","cov2","%ID","Description",""};

	Object[][] data;
	JTable table;

	String oldName1;
	String oldName2;

	String algorithmName;
	StructureAlignment algorithm;
	
	boolean isCE = true;
	UserConfiguration config;
	AtomCache cache ;

	String userPath ;
	String userChain;

	

	public static void main(String[] args){

		String file = "/tmp/results_4hhb.A.out";

		DBResultTable table = new DBResultTable();
		UserConfiguration config = WebStartMain.getDefaultConfig();
		config.setPdbFilePath("/Users/ap3/WORK/PDB/");
		config.setAutoFetch(true);
		config.setSplit(true);
		table.show(new File(file),config);
	}

	public DBResultTable(){
		oldName1 = "";
		oldName2 = "";
		userPath = null;
		userChain = null;
	}

	public void show(BufferedReader in, UserConfiguration config) throws IOException{
		String str;
		List<String[]> tmpdat = new ArrayList<String[]>();
		while ((str = in.readLine()) != null) {
			if ( str.startsWith("#")) {
				if ( str.startsWith("# algorithm:")) {
					String[] spl = str.split(":");
					if ( spl.length == 2) {
						algorithmName = spl[1];
						if (algorithmName.startsWith("jCE"))
							isCE = true;
						else
							isCE = false;
					}
					initAlgorithm(algorithmName);

				}

				else if ( str.startsWith("#param:file1=")){
					String path = str.substring(13);
					userPath = path.trim();
				}

				else if ( str.startsWith("#param:chain1=")){
					String chain = str.substring(14);
					userChain = chain.trim();
				} 
				
				else if ( str.startsWith("#param:scoring=")){
					try {
						String[] spl = str.split("=");
						String scoreS = spl[1];
						int scoring = Integer.parseInt(scoreS);
						if (algorithm != null){
							// scoring is a parameter of CE...
							ConfigStrucAligParams params = algorithm.getParameters();
							if ( params instanceof CeParameters){
								CeParameters ceParams = (CeParameters) params;
								ceParams.setScoringStrategy(scoring);
							}
						}
					} catch (Exception e){
						System.err.println("Unknown parameter can't read parameters from line: " + str);
						e.printStackTrace();
					}

				}
				continue;
			}
			String[] spl = str.split("\t");
			if ( spl.length != ceColumnNames.length -1) {
				System.err.println("wrong table width! " + spl.length + " should be: " + (ceColumnNames.length -1 ));
				System.err.println(str);
				continue;
			}
			tmpdat.add(spl);

		}
		in.close();
		
		Object[][] d = new Object[tmpdat.size()][ceColumnNames.length + 1];

		int i = -1; 
		for (String[] spl : tmpdat){
			
			i++;
			Object[] o = new Object[spl.length + 1];
			for ( int j=0; j< spl.length;j++){

				if (( j >= 2 && j <= 4)|| (j==9)) {
					o[j] = Double.parseDouble(spl[j]);
				}  else if (  j >4 && j< 10) {
				
					o[j] = Integer.parseInt(spl[j]);										    	    
				} else {
					o[j] = spl[j];
				}
			}

			o[spl.length ] = "Align";

			d[i] = o;

		}
		data = d;
		String[] columnNames = ceColumnNames;
		if ( ! isCE)
			columnNames = fatColumnNames;
		table = new JTable(data, columnNames);

		TableRowSorter<TableModel> sorter = new MyTableRowSorter(table.getModel());
		table.setRowSorter(sorter);
		//table.setAutoCreateRowSorter(true);

		JScrollPane scrollPane = new JScrollPane(table);
		table.setFillsViewportHeight(true);

		// take care of selections:
		table.setSelectionMode( ListSelectionModel.SINGLE_INTERVAL_SELECTION);
		table.getSelectionModel().addListSelectionListener(new RowListener());


		JFrame f = new JFrame();
		f.getContentPane().add(scrollPane);
		f.pack();
		f.setVisible(true);
		
	}
	
	public void show(File file, UserConfiguration config){
		this.config = config;

		cache = new AtomCache(config);
		try {
			BufferedReader in = new BufferedReader(new FileReader(file));
			show(in, config);
			
		} catch (IOException e) {
			e.printStackTrace();
		}

	}

	public void show(URL url, UserConfiguration config){
		this.config = config;

		cache = new AtomCache(config);
		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(url.openStream()));
			show(in, config);
			
		} catch (IOException e) {
			e.printStackTrace();
		}

	}

	
	private void initAlgorithm(String algorithmName) {
		try {
			algorithm = StructureAlignmentFactory.getAlgorithm(algorithmName);
		} catch (Exception e){
			e.printStackTrace();
			System.err.println("Can't guess algorithm from output. Using jCE as default...");
			try {
				algorithm = StructureAlignmentFactory.getAlgorithm(CeMain.algorithmName);
			} catch (Exception ex){
				ex.printStackTrace();
				return;
			}
		}

	}

	private void outputSelection() {
		StringBuffer output = new StringBuffer();
		output.append(String.format("Lead: %d, %d. ",
				table.getSelectionModel().getLeadSelectionIndex(),
				table.getColumnModel().getSelectionModel().
				getLeadSelectionIndex()));
		output.append("Rows:");
		for (int c : table.getSelectedRows()) {
			output.append(String.format(" %d", c));
		}

		output.append(". Columns:");
		for (int c : table.getSelectedColumns()) {
			output.append(String.format(" %d", c));
		}

		System.out.println(output.toString());
	}

	private class RowListener implements ListSelectionListener {
		public void valueChanged(ListSelectionEvent event) {
			if (event.getValueIsAdjusting()) {
				return;
			}            
			int row = table.getSelectionModel().getLeadSelectionIndex();
			String name1 = (String)table.getValueAt(row, 0);
			String name2 = (String)table.getValueAt(row, 1);

			if ( name1.equals(oldName1) && oldName2.equals(name2)){
				return;
			}
			System.out.println("recreating alignment of: " + name1 + " " + name2 + " using " + algorithmName);
			outputSelection();
			showAlignment(name1,name2);
			oldName1 = name1;
			oldName2 = name2;


		}
	}

	private void showAlignment( String name1, String name2){


		if ( algorithm == null) {
			initAlgorithm(null);
		}

		try {
			Structure structure1 = null;
			if ( name1.equals("CUSTOM")) {
				// user uploaded a custom PDB file...
				structure1 = loadCustomStructure(userPath,userChain);
			} else {						
				structure1 = cache.getStructure(name1);
			}
			Structure structure2 = cache.getStructure(name2);

			Atom[] ca1;
			Atom[] ca2;

			ca1 = StructureTools.getAtomCAArray(structure1);
			ca2 = StructureTools.getAtomCAArray(structure2);

			AFPChain afpChain;

			afpChain = algorithm.align(ca1, ca2);
			afpChain.setName1(name1);
			afpChain.setName2(name2);



			StructureAlignmentJmol jmol = StructureAlignmentDisplay.display(afpChain,ca1,ca2);

			//String result = afpChain.toFatcat(ca1, ca2);

			//String rot = afpChain.toRotMat();

			DisplayAFP.showAlignmentImage(afpChain, ca1,ca2,jmol);


		} catch (Exception e){
			e.printStackTrace();
		}
	}

	private Structure loadCustomStructure(String userPath2, String userChain2) throws StructureException{
		StructureIOFile reader = new PDBFileReader();
		Structure s = null;
		try {
			s = reader.getStructure(userPath2);
		} catch (IOException  e){

			//e.printStackTrace();
			throw new StructureException(e);
		}


		return StructureTools.getReducedStructure(s, userChain2);
	}

	public void actionPerformed(ActionEvent e) {



	}


}
