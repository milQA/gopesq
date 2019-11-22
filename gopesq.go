package gopesq

/*
#cgo LDFLAGS: -lm
#cgo CFLAGS: -I${SRCDIR}/include
#cgo CFLAGS: -I /usr/include
#include <cpesq/include/pesqmain.h>
#include <stdlib.h>


*/
import "C"
import (
	"unsafe"
)

func Match(file1, file2 string) float64 {
	const arrayLen = 4
	theCArray := [arrayLen]*C.char{
		C.CString(""),
		C.CString(file1),
		C.CString(file2),
		C.CString(""),
	}
	defer func() {
		for _, val := range theCArray {
			C.free(unsafe.Pointer(val))
		}
	}()

	return (float64(C.maiin(C.int(3), &theCArray[0])))
}
