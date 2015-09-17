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
 * Created on Jun 30, 2010
 * Author: ap3 
 *
 */
package org.biojava.nbio.structure.gui.util;

import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.align.util.UserConfiguration;
import org.biojava.nbio.structure.align.webstart.WebStartMain;

import java.util.ArrayList;
import java.util.List;

import javax.swing.*;

import java.awt.*;

/**
 * A Text Panel that allows the user to specify multiple structure 
 * identifiers, space separated.
 * 
 * @author Aleix Lafita
 * @since 4.1.1
 *
 */
public class SelectMultiplePanel extends JPanel {

	private static final long serialVersionUID = 757947454156959178L;

	JTextField input;

	public SelectMultiplePanel(){
		this(true);
	}

	public SelectMultiplePanel(boolean show2boxes){

		Box vBox = Box.createVerticalBox();

		input = new JTextField("1mbc 1hlb 1dlw 1ith.A 1thb.A 1kr7.A_0-110");
		Box b = getDomainPanel(input);
		vBox.add(b);

		this.add(vBox);
	}

	private Box getDomainPanel(JTextField f){

		JLabel l01 = new JLabel("Input structures:");

		Box hBox = Box.createHorizontalBox();
		hBox.add(Box.createGlue());
		hBox.add(l01);

		f.setMaximumSize(new Dimension(Short.MAX_VALUE,30));
		f.setToolTipText("Provide structure identifiers space separated.");

		hBox.add(Box.createVerticalGlue());
		hBox.add(f, BorderLayout.CENTER);
		hBox.add(Box.createGlue());

		return hBox;
	}

	public List<Structure> getStructures() throws StructureException {

		List<Structure> structures = new ArrayList<Structure>();

		for (String name:getNames()){
			structures.add(getStructure(name));
		}
		return structures;
	}

	public List<String> getNames() {

		List<String> names = new ArrayList<String>();

		String raw = input.getText().trim();
		String[] split = raw.split(" ");
		for (String name:split){
			if (name != null || name!="")
				names.add(name.trim());
		}
		return names;
	}

	private Structure getStructure(String name) throws StructureException{

		UserConfiguration config = WebStartMain.getWebStartConfig();

		AtomCache cache = new AtomCache(config);
		cache.setStrictSCOP(false);

		Structure s = null;
		try {
			s =	cache.getStructure(name);
			s.setName(name);
		} catch (Exception e){
			e.printStackTrace();
		}
		return s;
	}
}
