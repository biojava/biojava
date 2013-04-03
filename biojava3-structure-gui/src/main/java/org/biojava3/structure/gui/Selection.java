package org.biojava3.structure.gui;

import java.util.List;

import org.biojava.bio.structure.align.gui.jmol.AtomInfo;

public interface Selection {

	
	public void clear();
	public List<AtomInfo> getSelection();
	public void setSelection(List<AtomInfo> selection);
	
}
