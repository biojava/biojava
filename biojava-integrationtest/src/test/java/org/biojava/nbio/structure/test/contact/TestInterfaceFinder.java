package org.biojava.nbio.structure.test.contact;

import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureIO;
import org.biojava.nbio.structure.contact.AtomContact;
import org.biojava.nbio.structure.contact.AtomContactSet;
import org.biojava.nbio.structure.contact.InterfaceFinder;
import org.biojava.nbio.structure.contact.Pair;
import org.biojava.nbio.structure.contact.StructureInterface;
import org.biojava.nbio.structure.contact.StructureInterfaceList;
import org.junit.Test;

import java.io.IOException;
import java.util.HashSet;
import java.util.Set;

import static org.junit.Assert.assertEquals;

public class TestInterfaceFinder {

    @Test
    public void testGetAllInterfaces() throws StructureException, IOException {
        Structure s = StructureIO.getStructure("3hbx");

        InterfaceFinder finder = new InterfaceFinder(s);

        StructureInterfaceList list = finder.getAllInterfaces();

        assertEquals(12, list.size());

        Set<Pair<String>> unique = new HashSet<>();

        for (StructureInterface interf : list) {
            System.out.println("Interface " + interf.getMoleculeIds());
            AtomContactSet set = interf.getContacts();
            System.out.println("Number of contacts: " + set.size());

            unique.add(interf.getMoleculeIds());

        }
        assertEquals(12, unique.size());
    }
}
