package org.biojava.nbio.core.util;

import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.time.LocalDate;
import java.time.Month;
import java.time.temporal.Temporal;
import java.util.ArrayList;
import java.util.List;

import org.junit.jupiter.api.Nested;
import org.junit.jupiter.api.Test;

class EqualsTest {

    @Test
    void equalsInt(){
        assertTrue(Equals.equal(1, 1));
        assertTrue(Equals.equal(Integer.MAX_VALUE, Integer.MAX_VALUE));
        assertTrue(Equals.equal(Integer.valueOf(1), Integer.valueOf(1)));
        assertFalse(Equals.equal(1, 2));
    }

    void equalsBool(){
        assertTrue(Equals.equal(true, true));
        assertTrue(Equals.equal(Boolean.valueOf(true), Boolean.valueOf(true)));
        assertFalse(Equals.equal(true, false));
    }

    void equalsLong(){
        assertTrue(Equals.equal(1L, 1L));
        assertTrue(Equals.equal(Long.valueOf(1L), Long.valueOf(1L)));
        assertTrue(Equals.equal(Long.MAX_VALUE, Long.MAX_VALUE));
        assertFalse(Equals.equal(1L, 1L));
    }

    @Nested
    class ObjectEquals {
        Object o1 = new Object();
        Object o1Ref = o1;
        Object o2 = new Object();

        @Test
        void twoNullsAreEqual(){
            assertTrue(Equals.equal(null, null));
        }
        @Test
        void objectWithNullIsNotEqual(){
            assertFalse(Equals.equal(o1, null));
            assertFalse(Equals.equal(null, o1));
        }

        @Test
        void identicalObjectIsEquals(){
            assertTrue(Equals.equal(o1, o1));
            assertTrue(Equals.equal(o1, o1Ref));
            assertFalse(Equals.equal(o1, o2));
        }

        @Test
        void equalsBasedOnProperties(){
            LocalDate date = LocalDate.of(2021, Month.APRIL, 21);
            LocalDate sameDate =  LocalDate.of(2021, Month.APRIL, 21);
            LocalDate differentDate =  LocalDate.of(2022, Month.APRIL, 21);
            assertTrue(Equals.equal(date, sameDate));
            assertFalse(Equals.equal(date, differentDate));
        }
    }

    @Nested
    class ClassEquals {
        LocalDate nowDate = LocalDate.now();
        LocalDate different = nowDate.plusDays(5);
        @Test
        void identicalClassesAreEqual() {  
            assertTrue(Equals.classEqual(nowDate, different));
        }
        @Test
        void classComparisonIsByActualTypeNotReferenceType() {
            Temporal temporal = (Temporal) nowDate;
            assertTrue(Equals.classEqual(temporal, different));
        }

        @Test
        void genericsAreIgnored() {
            List<String> listOfStrings = new ArrayList<>();
            List<Integer> listOfInts = new ArrayList<>();
            assertTrue(Equals.classEqual(listOfStrings, listOfInts));
        }

        class ASuperclass {
            Integer a, b = 0;
        }
        class ASubclass extends ASuperclass {
            Integer c, d = 0;
        }
        @Test
        void membersOfClassHierarchyAreNotEqual() {
            ASuperclass superObject = new ASuperclass();
            ASuperclass subObject = new ASubclass();
            assertFalse(Equals.classEqual(superObject, subObject));
            assertFalse(Equals.classEqual(subObject, superObject));
        }   
    }  
}
