package org.biojava.nbio.core.util;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;

import org.junit.jupiter.api.Test;

class SoftHashMapTest {

    static class TestObject  {
        private String name;
        public TestObject(String string) {
            this.name=string;
        }
        String getName(){
            return name;
        }
    }
    
    @Test
    void basicMapOperations() throws InterruptedException{
        
        SoftHashMap<String, TestObject> map = new SoftHashMap<>(1);
        TestObject s1= new TestObject("1");
        TestObject s2= new TestObject("2");
        TestObject s3= new TestObject("3");
    
        map.put("1", s1);
        map.put("2", s2);
        map.put("3", s3);
        assertEquals(3, map.size());

        map.put("3", new TestObject("4"));
        assertEquals(3, map.size());

        assertEquals(s1, map.remove("1"));
        assertEquals(2, map.size());

        map.clear();
        assertEquals(0, map.size());
    }
    @Test
    void manyMapOperationsAreUnsupported() throws Exception{
        SoftHashMap<String, TestObject> map = new SoftHashMap<>(1);
        TestObject s1= new TestObject("1");
        map.put("1", null);
        // these all use entrySet internally and throw USOException
        assertThrows(UnsupportedOperationException.class, ()->map.containsValue(s1));
        assertThrows(UnsupportedOperationException.class, ()->map.containsKey("1"));
        assertThrows(UnsupportedOperationException.class, ()->map.values().iterator());
        assertThrows(UnsupportedOperationException.class, ()->map.getOrDefault("1", new TestObject("2")));     
    }
        
}
