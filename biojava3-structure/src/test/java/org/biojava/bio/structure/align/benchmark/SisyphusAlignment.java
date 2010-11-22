/**
 * 
 */
package org.biojava.bio.structure.align.benchmark;

//import java.io.IOException;
//import java.util.ArrayList;
//import java.util.HashMap;
//import java.util.Iterator;
//import java.util.List;
//import java.util.Map;
//import java.util.Set;
//
//import org.biojava.bio.Annotation;
//import org.biojava.bio.program.das.dasalignment.AFPChainDASAlignmentConverter;
//import org.biojava.bio.program.das.dasalignment.Alignment;
//import org.biojava.bio.structure.Chain;
//import org.biojava.bio.structure.ChainImpl;
//import org.biojava.bio.structure.Group;
//import org.biojava.bio.structure.Structure;
//import org.biojava.bio.structure.StructureException;
//import org.biojava.bio.structure.StructureImpl;
//import org.biojava.bio.structure.align.gui.StructureAlignmentDisplay;
//import org.biojava.bio.structure.align.gui.jmol.StructureAlignmentJmol;
//import org.biojava.bio.structure.align.model.AFPChain;
//import org.biojava.bio.structure.align.util.AtomCache;
//import org.biojava.bio.structure.io.PDBFileReader;
//import org.biojava.bio.structure.jama.Matrix;
//import org.biojava.dasobert.das.AlignmentParameters;
//import org.biojava.dasobert.das.AlignmentThread;
//import org.biojava.dasobert.dasregistry.Das1Source;
//import org.biojava.dasobert.eventmodel.AlignmentEvent;
//import org.biojava.dasobert.eventmodel.AlignmentListener;
//import org.biojava3.structure.gui.JmolViewerImpl;
//import org.jmol.api.JmolViewer;

/**
 * @author Spencer Bliven <sbliven@ucsd.edu>
 *
 */
public class SisyphusAlignment {//implements MultipleAlignmentParser {
//	private Das1Source dasSource;
//	
//	public SisyphusAlignment() {
//		dasSource = new Das1Source();
//
//		dasSource.setUrl("http://sisyphus.mrc-cpe.cam.ac.uk/sisyphus/das/alignments/");
//	}
//
//	public AFPChain getAlignment(String sisyID, Structure s1, Structure s2 ) {
//		//String pdb1="/Users/blivens/pdb/1g7s.pdb";
//		//String pdb2="/Users/blivens/pdb/1kmh.pdb";
//		//sisyFrag1=809-782;
//		//sisyFrag2=794-782;
//
//		//reference1 = loadStructure(pdb1);
//		//reference2 = loadStructure(pdb2);
//
//		AlignmentParameters params = new AlignmentParameters();
//		params.setDasSource(dasSource);
//		params.setQuery(sisyID);
//
//		AlignmentThread t = new AlignmentThread(params);
//		MyAlignmentListener al = new MyAlignmentListener();
//		t.addAlignmentListener(al);
//		t.run();
//		AFPChain afpChain = null;
//		synchronized (t) {
//			try {
//				while(!al.isFinished()) {
//					t.wait();
//				}
//				Alignment[] alignments = al.getAlignments();
//				for( Alignment align:alignments) {
//					AFPChain afp = AFPChainDASAlignmentConverter.convertAlignmentToAFPChain(align, s1, s2);
//					if(afp != null ) {
//						afpChain = afp;
//					}
//				}
//				
//			} catch (InterruptedException e) {
//				e.printStackTrace();
//				return null;
//			}
//		}
//		
//		return afpChain;
//	}
//	/**
//	 * Creates a Structure object from a PDB file
//	 * @param pathToPDBFile
//	 * @return
//	 *
//	public static Structure loadStructure(String pathToPDBFile){
//		PDBFileReader pdbreader = new PDBFileReader();
//
//		Structure structure = null;
//		try{
//			structure = pdbreader.getStructure(pathToPDBFile);
//			System.out.println(structure);
//		} catch (IOException e) {
//			e.printStackTrace();
//		}
//		return structure;
//	}*/
// 
//	private class MyAlignmentListener implements AlignmentListener {
//		private boolean parsingFinished;
//		private Alignment[] alignments;
//		
//		public MyAlignmentListener() {
//			parsingFinished=false;
//			alignments = null;
//		}
//		
//		public boolean isFinished() {
//			return parsingFinished;
//		}
//		
//		public Alignment[] getAlignments() {
//			return alignments;
//		}
//		
//		@Override
//		public void noAlignmentFound(AlignmentEvent e)
//		{
//			parsingFinished=true;
//			alignments = e.getAllAlignments();
//		}
//
//		@Override
//		public void newAlignment(AlignmentEvent e)
//		{
//			parsingFinished = true;
//			alignments = e.getAllAlignments();
//
//			
//			/*try {
//				System.out.println("Filtering");
//				reference1 = getStructureRanges(sisyFrag1,reference1,a);
//				reference2 = getStructureRanges(sisyFrag2,reference2,a);
//
//				//TODO StructureAlignmentJmol jmol = StructureAlignmentDisplay.display(afpChain,ca1,ca2,hetatms,nucs, hetatms2, nucs2);
//
//				System.out.println("Reference 1");
//				for(Chain c : reference1.getChains()) {
//					for(Group g :c.getGroups()) {
//						System.out.format("%s:%s\t%s\n", g.getPDBCode(), c.getName(), g.getPDBName() );
//					}
//				}
//				System.out.println("Reference 2");
//				for(Chain c : reference2.getChains()) {
//					for(Group g :c.getGroups()) {
//						System.out.format("%s:%s\t%s\n", g.getPDBCode(), c.getName(), g.getPDBName() );
//					}
//				}
//			} catch (StructureException e1) {
//				// TODO Auto-generated catch block
//				e1.printStackTrace();
//			}*/
//		}
//
//		@Override
//		public void clearAlignment()
//		{
//			// This seems to never get called.
//			throw new UnsupportedOperationException("clearAlignment Handler not implemented");
//		}
//	}
//
//	/** return the structure object restricted to the ranges given in the ALignment
//	 * if no ranges give, returns null
//	 * @param pos
//	 * @param s
//	 * @return null; if no ranges given, otherwise restricts the
//		structure to the ranges.
//	 * @throws StructureException
//	 *
//	private static Structure getStructureRanges(int pos , Structure s, Alignment alignment)
//	throws StructureException{
//
//		if (alignment == null)
//			return null;
//
//		// check if the alignment has an object detail "region" for this
//		// if yes = restrict the used structure...
//		Annotation[] objects = alignment.getObjects();
//		Annotation object = objects[pos];
//		//System.out.println(object);
//		
//		// Keeps track of residues to keep around. Map ChainID->[PDB Codes]
//		Map<String,List<String>> protectionMap = new HashMap<String,List<String>>();
//
//		// for each alignObject
//		if ( object.containsProperty("details")){
//
//			List<Annotation> details = (List<Annotation>) object.getProperty("details");
//
//			// for each alignObjectDetail
//			for ( int det = 0 ; det< details.size();det++) {
//				Annotation detanno = details.get(det);
//				String property = (String)detanno.getProperty("property");
//				if (! property.equals("region"))
//					continue; //ignore none-region details
//
//				//get the region string
//				//Regions are formatted like A:1-999
//				String detail = (String) detanno.getProperty("detail");
//				//System.out.println("got detail " + detail);
//
//
//				// split up the structure and add the region to the new structure...
//				int cpos = detail.indexOf(":");
//				String chainId = " ";
//
//				if ( cpos > 0) {
//					chainId = detail.substring(0,cpos);
//					detail  = detail.substring(cpos+1,detail.length());
//				} else {
//					detail = detail.substring(1,detail.length());
//				}
//
//				//System.out.println(detail + " " + cpos + " " + chainId);
//
//				String[] spl = detail.split("-");
//
//				if ( spl.length != 2)
//					continue;
//				String start = spl[0];
//				String end   = spl[1];
//				//System.out.println("start " + start + " end " + end);
//
//
//				List<String> prot = protectionMap.get(chainId);
//
//				List<String> protectedResidues;
//				if ( prot == null)
//					protectedResidues = new ArrayList<String>();
//				else {
//					protectedResidues = prot;
//				}
//
//				Chain c = s.getChainByPDB(chainId);
//
//				//TODO: Q: do we need to do  all groups or only the aminos??
//				List<Group> groups = c.getGroups();
//
//				boolean known =false;
//
//				for(Group g: groups){
//					
//					if (g.getPDBCode().equals(start)){
//						known = true;
//					}
//
//					Group n = (Group) g.clone();
//
//					String code = n.getPDBCode();
//
//					if ( known)
//						protectedResidues.add(code);
//
//					if (g.getPDBCode().equals(end)){
//						known = false;
//					}
//
//				}
//				protectionMap.put(chainId,protectedResidues);
//
////				/*Group[] groups = c.getGroupsByPDB(start,end);
////
////		                Chain nc = new ChainImpl ();
////		                nc.setName(chainId);
////		                boolean knownChain = false;
////		                try {
////		                nc = newStruc.findChain(chainId);
////		                knownChain = true;
////
////		                } catch (Exception e){}
////
////		                for (int g=0;g<groups.length;g++){
////		                Group gr = groups[g];
////		                nc.addGroup(gr);
////
////
////		                }
////				 *
//			}
//		}
//
//		System.out.println(protectionMap.toString());
//		
//		// o.k. here we got all residues in the protectedMap
//		// now we build up the new structure ...
//
//		Structure newStruc = new StructureImpl();
//		newStruc.setPDBCode(s.getPDBCode());
//		newStruc.setPDBHeader(s.getPDBHeader());
//
//		for( String chainId : protectionMap.keySet() ) {
//			List<String> residues = protectionMap.get(chainId);
//
//			Chain origChain = s.findChain(chainId);
//			Chain newChain  = new ChainImpl();
//			newChain.setName(chainId);
//
//			//TODO: Q: do we need to do  all groups or only the aminos??
//			List<Group> origGroups = origChain.getGroups();
//			for(Group orig : origGroups) {
//				Group n = (Group) orig.clone();
//				if ( ! residues.contains(orig.getPDBCode())){
//					// this group has been
//					n.clearAtoms();
//
//				}
//				newChain.addGroup(n);
//			}
//			newStruc.addChain(newChain);
//		}
//
//		if (newStruc.size() > 0){
//			return newStruc;
//		} else
//			return null;
//
//	}
//	*/
//
//	/**
//	 * @param args
//	 */
//	public static void main(String[] args) {
//		String sisyID="AL10050412";
//		String pdb1="1g7s.A";
//		String pdb2="1kmh.A";
//		//sisyFrag1=809-782;
//		//sisyFrag2=794-782;
//		AtomCache cache = new AtomCache("/tmp/", true);
//
//		Structure s1, s2;
//		try {
//			s1 = cache.getStructure(pdb1);
//			s2=cache.getStructure(pdb2);
//		} catch (StructureException e) {
//			e.printStackTrace();
//			return;
//		} catch (IOException e) {
//			e.printStackTrace();
//			return;
//		}
//
//		SisyphusAlignment sisy = new SisyphusAlignment();
//		AFPChain alignment = sisy.getAlignment(sisyID,s1,s2);
// 		System.out.println("Got an AFPChains");
//
//	}
//
//
}
