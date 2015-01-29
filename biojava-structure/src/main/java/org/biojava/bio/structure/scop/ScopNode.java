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

package org.biojava.bio.structure.scop;

import java.io.Serializable;
import java.util.List;

import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;
import javax.xml.bind.annotation.XmlRootElement;

@XmlRootElement(name = "ScopNode", namespace ="http://source.rcsb.org")
@XmlAccessorType(XmlAccessType.PUBLIC_MEMBER)
public class ScopNode implements Serializable
{

	/**
	 * 
	 */
	private static final long serialVersionUID = 1187083944488580995L;

	int sunid;
	int parentSunid;
	List<Integer> children;

	public ScopNode(){

	}



	@Override
	public String toString()
	{
		return "ScopNode [children=" + children + ", parentSunid=" + parentSunid + ", sunid=" + sunid + "]";
	}



	public int getSunid()
	{
		return sunid;
	}
	public void setSunid(int sunid)
	{
		this.sunid = sunid;
	}
	public int getParentSunid()
	{
		return parentSunid;
	}
	public void setParentSunid(int parentSunid)
	{
		this.parentSunid = parentSunid;
	}
	public List<Integer> getChildren()
	{
		return children;
	}
	public void setChildren(List<Integer> children)
	{
		this.children = children;
	}


}
