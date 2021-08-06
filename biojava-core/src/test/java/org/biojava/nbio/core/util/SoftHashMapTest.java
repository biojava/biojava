package org.biojava.nbio.core.util;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.Map;

import org.junit.jupiter.api.Disabled;
import org.junit.jupiter.api.Test;
/** 
 * Includes a disabled test that asserts behaviour of collecting
 * SoftReferences, run using -Xmx=5M to expose this behaviour.
 */
class SoftHashMapTest {

    static class TestObject  {
        /*
        *Create an object occupying negligible memory
        */
        static TestObject small(String name){
            return new TestObject(name, 100);
        }

        /*
        *Create a test object occupying significant memory(100kB)
        */
        static TestObject large(String name){
            return new TestObject(name, 100_000);
        }
        private String name;
        private int [] internalArray = null; 
        public TestObject(String string, int capacity) {
            this.name=string;
            this.internalArray = new int [capacity];
        }
        String getName(){
            return name;
        }
        public String toString(){
            return name;
        }
    }

    // This test needs to be run with restricted memory in order
    // to expose the behaviour of SoftHashMap in deleting entries 
    // when under memory pressure. By setting -Xmx=5M the test can
    // assert that entries are deleted. We don't want to risk throwing
    // OOM errors during normal test execution so this is disabled
    @Test
    @Disabled("requires to run in conditions nearly throwing an OOM")
    void softMapRemovesRefsToSaveMemory() throws InterruptedException{
       
        // Using a regular Map with hard references will probably
        // cause an OOM error if running with -Xmx=5M. Uncomment this
        // and comment out the next line to observe this.
        // Map<String, TestObject> map =new HashMap<>(1);

        // set the maximum number of hard references to 1 (minimum)
        // to expose behaviour of soft references better. 
        Map<String, TestObject> map = new SoftHashMap<>(1);
        int totalPuts =5;
        for (int i = 0; i < totalPuts; i++) {

            TestObject myObject = TestObject.large(""+i);
            map.put(myObject.getName(),myObject);
            //allocate a little slowly
            // enables GC time to work
            Thread.sleep(10);
        }
        int nonNullValues = countNonNullMapReferences(map, totalPuts);
        // some but not all references should be removed.
        assertTrue(nonNullValues > 0 && nonNullValues < totalPuts);
    }

    private int countNonNullMapReferences(Map<String, TestObject> map, int totalPuts) {
        try {
            //sleep a little in case if finalizers are currently running
            Thread.sleep(1000);
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
        // we can't iterate over map as keySet() isn't implemented
        int nonNullValues = 0;
        for (int i = 0; i< totalPuts; i++) {
            if(map.get("" + i ) != null) {
                nonNullValues++;
            }
        } 
        return nonNullValues;
    }
    
    @Test
    void basicMapOperations() throws InterruptedException{
        
        SoftHashMap<String, TestObject> map = new SoftHashMap<>(1);
        TestObject s1= TestObject.small("1");
        TestObject s2= TestObject.small("2");
        TestObject s3= TestObject.small("3");
    
        map.put("1", s1);
        map.put("2", s2);
        map.put("3", s3);
        assertEquals(3, map.size());

        map.put("3", TestObject.small("4"));
        assertEquals(3, map.size());

        assertEquals(s1, map.remove("1"));
        assertEquals(2, map.size());

        map.clear();
        assertEquals(0, map.size());
    }
    @Test
    void manyMapOperationsAreUnsupported() throws Exception{
        SoftHashMap<String, TestObject> map = new SoftHashMap<>(1);
        TestObject s1= TestObject.small("1");
        map.put("1", null);
        // these all use entrySet internally and throw USOException
        assertThrows(UnsupportedOperationException.class, ()->map.containsValue(s1));
        assertThrows(UnsupportedOperationException.class, ()->map.containsKey("1"));
        assertThrows(UnsupportedOperationException.class, ()->map.values().iterator());
        assertThrows(UnsupportedOperationException.class, ()->map.getOrDefault("1", TestObject.small("2")));     
    }
        
}
