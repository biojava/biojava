package demo;



import java.util.List;

import javax.vecmath.AxisAngle4d;
import javax.vecmath.Vector3d;

import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.contact.AtomContact;
import org.biojava.bio.structure.contact.Pair;
import org.biojava.bio.structure.contact.StructureInterface;
import org.biojava.bio.structure.contact.StructureInterfaceList;
import org.biojava.bio.structure.xtal.CrystalBuilder;
import org.biojava.bio.structure.xtal.CrystalTransform;
import org.biojava.bio.structure.xtal.SpaceGroup;
import org.biojava3.structure.StructureIO;



public class DemoCrystalInterfaces {

	
	private static final double BSATOASA_CUTOFF = 0.95;
	private static final double MIN_ASA_FOR_SURFACE = 5;
	private static final int CONSIDER_COFACTORS = 40; // minimum number of atoms for a cofactor to be considered, if -1 all ignored
	
	
	private static final double CUTOFF = 5.5; 
	
	private static final int N_SPHERE_POINTS = 3000;
	
	private static final double MIN_AREA_TO_KEEP = 35;
	
	private static final int NTHREADS = Runtime.getRuntime().availableProcessors();
	
	private static final boolean DEBUG = true;
	
	private static final double CLASH_DISTANCE = 1.5; 
	
	
	/**
	 * @param args
	 */
	public static void main(String[] args) throws Exception {

		
		String pdbCode = "1smt";
		
		AtomCache cache = new AtomCache();
		cache.setUseMmCif(true);
		
		StructureIO.setAtomCache(cache); 
		
		Structure structure = StructureIO.getStructure(pdbCode);
			

		System.out.println(structure.getPDBCode());
		
		
		SpaceGroup sg = structure.getCrystallographicInfo().getSpaceGroup();
		
		System.out.println(sg.getShortSymbol()+" ("+sg.getId()+")");
		System.out.println("Symmetry operators: "+sg.getNumOperators());
		
		System.out.println("Calculating possible interfaces... (using "+NTHREADS+" CPUs for ASA calculation)");
		long start = System.currentTimeMillis();
		
		CrystalBuilder cb = new CrystalBuilder(structure);
		cb.setDebug(DEBUG); 
		
		
		StructureInterfaceList interfaces = cb.getUniqueInterfaces(CUTOFF);
		interfaces.setDebug(DEBUG);
		interfaces.calcAsas(N_SPHERE_POINTS, NTHREADS, CONSIDER_COFACTORS);
		interfaces.removeInterfacesBelowArea(MIN_AREA_TO_KEEP);

		
		//interfaces.initialiseClusters(pdb, CLUSTERING_CUTOFF, MINATOMS_CLUSTERING, "CA");
		
		long end = System.currentTimeMillis();
		long total = (end-start)/1000;
		System.out.println("Total time for interface calculation: "+total+"s");
		
		System.out.println("Total number of interfaces found: "+interfaces.size());

		for (int i=0;i<interfaces.size();i++) {
			StructureInterface interf = interfaces.get(i+1);			
			
			String infiniteStr = "";
			if (interf.isInfinite()) infiniteStr = " -- INFINITE interface";
			System.out.println("\n##Interface "+(i+1)+" "+
					interf.getCrystalIds().getFirst()+"-"+
					interf.getCrystalIds().getSecond()+infiniteStr);
			// warning if more than 10 clashes found at interface
			List<AtomContact> clashing = interf.getContacts().getContactsWithinDistance(CLASH_DISTANCE);
			if (clashing.size()>10) 
				System.out.println(clashing.size()+" CLASHES!!!");
			
			CrystalTransform transf1 = interf.getTransforms().getFirst();
			CrystalTransform transf2 = interf.getTransforms().getSecond();
			
			System.out.println("Transf1: "+SpaceGroup.getAlgebraicFromMatrix(transf1.getMatTransform())+
					". Transf2: "+SpaceGroup.getAlgebraicFromMatrix(transf2.getMatTransform()));
	 		
			int foldType = sg.getAxisFoldType(transf2.getTransformId());
			AxisAngle4d axisAngle = sg.getRotAxisAngle(transf2.getTransformId());
			
			String screwStr = "";
			if (transf2.getTransformType().isScrew()) {
				Vector3d screwTransl = 
						transf2.getTranslScrewComponent();
				screwStr = " -- "+transf2.getTransformType().getShortName()+" with translation "+
				String.format("(%5.2f,%5.2f,%5.2f)",screwTransl.x,screwTransl.y,screwTransl.z);

			}
			
			
			System.out.println(" "+foldType+"-fold on axis "+String.format("(%5.2f,%5.2f,%5.2f)",axisAngle.x,axisAngle.y,axisAngle.z)+screwStr);
			
			System.out.println("Number of contacts: "+interf.getContacts().size());
			//System.out.println("Number of contacting atoms (from both molecules): "+interf.getNumAtomsInContact());
			Pair<List<Group>> cores = interf.getCoreResidues(BSATOASA_CUTOFF, MIN_ASA_FOR_SURFACE);
			System.out.println("Number of core residues at "+String.format("%4.2f", BSATOASA_CUTOFF)+
					" bsa to asa cutoff: "+
					cores.getFirst().size()+" "+
					cores.getSecond().size());
			System.out.printf("Interface area: %8.2f\n",interf.getTotalArea());
			
		}
		


	}

}
