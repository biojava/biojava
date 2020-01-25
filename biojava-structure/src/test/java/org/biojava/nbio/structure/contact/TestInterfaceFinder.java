package org.biojava.nbio.structure.contact;

import org.biojava.nbio.structure.AminoAcidImpl;
import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.AtomImpl;
import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.ChainImpl;
import org.biojava.nbio.structure.EntityInfo;
import org.biojava.nbio.structure.Group;
import org.biojava.nbio.structure.ResidueNumber;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureImpl;
import org.biojava.nbio.structure.StructureTools;
import org.junit.Test;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import static org.junit.Assert.*;

import javax.vecmath.Point3d;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class TestInterfaceFinder {

    private static final Logger logger = LoggerFactory.getLogger(TestInterfaceFinder.class);

    @Test
    public void testGetAllInterfaces() {
        Structure s = mockStructure();
        InterfaceFinder finder = new InterfaceFinder(s);

        StructureInterfaceList list = finder.getAllInterfaces();

        assertEquals(3, list.size());

        Set<Pair<String>> unique = new HashSet<>();

        for (StructureInterface interf : list) {
            System.out.println("Interface " + interf.getMoleculeIds());
            AtomContactSet set = interf.getContacts();
            for (AtomContact c : set)
                System.out.println(c.getPair() +" - " + c.getDistance());

            unique.add(interf.getMoleculeIds());

        }
        assertEquals(3, unique.size());
    }

    /**
     * Create a mock structure with 2 entities 1 (chains A, B) and 2 (chain C).
     * @return a structure
     */
    private Structure mockStructure() {
        Structure structure = new StructureImpl();
        EntityInfo entity1 = new EntityInfo();
        entity1.setMolId(1);
        EntityInfo entity2 = new EntityInfo();
        entity2.setMolId(2);
        structure.addEntityInfo(entity1);
        structure.addEntityInfo(entity2);

        Chain chainA = new ChainImpl();
        chainA.setId("A");
        chainA.setName("A");
        Chain chainB = new ChainImpl();
        chainB.setId("B");
        chainB.setName("B");
        entity1.addChain(chainA);
        entity1.addChain(chainB);
        Chain chainC = new ChainImpl();
        chainC.setId("C");
        chainC.setName("C");
        entity2.addChain(chainC);

        structure.addChain(chainA);
        structure.addChain(chainB);
        structure.addChain(chainC);

        // entity 1: chain A 10 observed residues, chain B 9 observed residues (first unobserved)
        List<Group> aGroups = getGroupList(10, "ALA", chainA, new Point3d(0,0,0));
        chainA.setAtomGroups(new ArrayList<>(aGroups));
        chainA.setSeqResGroups(aGroups);
        chainA.setEntityInfo(entity1);

        List<Group> bGroups = getGroupList(10, "ALA", chainB, new Point3d(4, 0, 0));
        chainB.setAtomGroups(new ArrayList<>(bGroups.subList(1,10)));
        chainB.setSeqResGroups(bGroups);
        chainB.setEntityInfo(entity1);

        List<Group> cGroups = getGroupList(20, "GLY", chainC, new Point3d(0, 4, 0));
        chainC.setAtomGroups(new ArrayList<>(cGroups));
        chainC.setSeqResGroups(cGroups);
        chainC.setEntityInfo(entity2);

        return structure;
    }

    private List<Group> getGroupList(int size, String type, Chain chain, Point3d center) {
        List<Group> list = new ArrayList<>();
        double offsetx = 0;
        double offsety = 0;
        double offsetz = 0;
        for (int i=0;i<size;i++) {
            Group g = new AminoAcidImpl();
            g.setPDBName(type);
            g.setResidueNumber(new ResidueNumber(chain.getId(), i+1, null));
            chain.addGroup(g);
            Atom a = new AtomImpl();
            a.setName(StructureTools.CA_ATOM_NAME);
            a.setX(center.x + offsetx);
            a.setY(center.y + offsety);
            a.setZ(center.z + offsetz);
            g.addAtom(a);
            list.add(g);

            if (i%3 == 0) offsetx += 1;
            if (i%3 == 1) offsety += 1;
            if (i%3 == 2) offsetz += 1;
        }

        return list;
    }

}
