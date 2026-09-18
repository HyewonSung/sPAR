	.file	"spqlios-fft-avx.s"
	.section	.text.unlikely,"ax",@progbits
.LCOLDB0:
	.text
.LHOTB0:
	.p2align 4,,15
	.globl	fft
	.type	fft, @function
fft:
.LFB0:
	.cfi_startproc







	/* Save registers */
	pushq       %r10
	pushq       %r11
	pushq       %r12
	pushq       %r13
	pushq       %r14
	pushq       %rbx
	
	/* Permute registers for better variable names */
	movq        %rdi, %rax
	movq        %rsi, %rdi      /* rdi: base of the real data CONSTANT */
	
	/* Load struct FftTables fields */
	movq         0(%rax), %rdx  /* rdx: n (logical Size of FFT  = a power of 2, must be at least 4) */
	movq         8(%rax), %r8   /* r8: Base address of trigonometric tables array (CONSTANT) */
	



	movq	%rdx, %r9
	shr	$2,%r9              /* r9: ns4 CONSTANT */
	leaq	(%rdi,%r9,8),%rsi   /* rsi: base of imaginary data CONSTANT */






























	vmovapd     size4negation0, %ymm15
	vmovapd     size4negation1, %ymm14
	vmovapd     size4negation2, %ymm13
	vmovapd     size4negation3, %ymm12
	
	movq	$0,%rax	/* rax: block */
	movq	%rdi,%r10
	movq	%rsi,%r11
fftsize2loop:
	vmovapd (%r10),%ymm0 /* r0 r1 r2 r3 */
	vmovapd (%r11),%ymm1 /* i0 i1 i2 i3 */
	vshufpd $0,%ymm0,%ymm0,%ymm2  /* r0 r0 r2 r2 */
	vshufpd $15,%ymm0,%ymm0,%ymm3 /* r1 r1 r3 r3 */
	vshufpd $0,%ymm1,%ymm1,%ymm4  /* i0 i0 i2 i2 */
	vshufpd $15,%ymm1,%ymm1,%ymm5 /* i1 i1 i3 i3 */
	vfmadd231pd %ymm3,%ymm12,%ymm2 /* (r0 r0 r2 r2) + (r1 -r1 r3 -r3) */
	vfmadd231pd %ymm5,%ymm12,%ymm4 /* (i0 i0 i2 i2) + (i1 -i1 i3 -i3) */
	vmovapd %ymm2,(%r10)
	vmovapd %ymm4,(%r11)
	/* end of loop */
        leaq 32(%r10),%r10
        leaq 32(%r11),%r11
        addq $4,%rax
	cmpq %r9,%rax
	jb fftsize2loop
	


































	movq	$0, %rax
	movq	%rdi,%r10
	movq	%rsi,%r11
fftsize4loop:
	vmovapd (%r10),%ymm0 /* r0 r1 r2 r3 */
	vmovapd (%r11),%ymm1 /* i0 i1 i2 i3 */
	vperm2f128 $32,%ymm0,%ymm0,%ymm4 /* r0 r1 r0 r1 */
	vperm2f128 $32,%ymm1,%ymm1,%ymm5 /* i0 i1 i0 i1 */
	vperm2f128 $49,%ymm0,%ymm0,%ymm6 /* r2 r3 r2 r3 */
	vperm2f128 $49,%ymm1,%ymm1,%ymm7 /* i2 i3 i2 i3 */
	vshufpd $10,%ymm7,%ymm6,%ymm8    /* r2 i3 r2 i3 */
	vshufpd $10,%ymm6,%ymm7,%ymm9    /* i2 r3 i2 r3 */
	vfmadd231pd %ymm8,%ymm13,%ymm4   /* (r0 r1 r0 r1) + (r2 i3 -r2 -i3) */
	vfmadd231pd %ymm9,%ymm14,%ymm5   /* (i0 i1 i0 i1) + (i2 -r3 -i2 r3) */
	vmovapd %ymm4,(%r10)
  	vmovapd %ymm5,(%r11)
        /* end of loop */
        leaq 32(%r10),%r10
        leaq 32(%r11),%r11
        addq $4,%rax
	cmpq %r9,%rax
	jb fftsize4loop


































	movq %r8,%rdx /* rdx: cur_tt */
	movq $4,%rax /* rax: halfnn */
ffthalfnnloop:
	movq $0,%rbx /* rbx: block */
fftblockloop:
	leaq (%rdi,%rbx,8),%r10 /* re0 pointer */
	leaq (%rsi,%rbx,8),%r11 /* im0 pointer */
	leaq (%r10,%rax,8),%r12 /* re1 pointer */
	leaq (%r11,%rax,8),%r13 /* im1 pointer */
	movq %rdx,%r14          /* tcs pointer */
	movq $0,%rcx /* rcx: off */
fftoffloop:
	vmovapd (%r10),%ymm0 /* re0 */
	vmovapd (%r11),%ymm1 /* im0 */
	vmovapd (%r12),%ymm2 /* re1 */
	vmovapd (%r13),%ymm3 /* im1 */
	vmovapd (%r14),%ymm4 /* cos */
	vmovapd 32(%r14),%ymm5 /* sin */
	vmulpd	%ymm2,%ymm4,%ymm6 /* re1.cos */
	vmulpd	%ymm2,%ymm5,%ymm7 /* re1.sin */
        vfnmadd231pd %ymm3,%ymm5,%ymm6 /* re2 = re1.cos - im1.sin */ 
        vfmadd231pd %ymm3,%ymm4,%ymm7  /* im2 = re1.sin + im1.cos */ 
	vsubpd	%ymm6,%ymm0,%ymm2 /* re0 - re2 */
	vsubpd	%ymm7,%ymm1,%ymm3 /* im0 - im2 */
	vaddpd	%ymm6,%ymm0,%ymm0 /* re0 + re2 */
	vaddpd	%ymm7,%ymm1,%ymm1 /* im0 + im2 */
	vmovapd %ymm0,(%r10)
	vmovapd %ymm1,(%r11)
	vmovapd %ymm2,(%r12)
	vmovapd %ymm3,(%r13)
        /* end of off loop */
    	leaq 	32(%r10),%r10
    	leaq	32(%r11),%r11
    	leaq 	32(%r12),%r12
    	leaq 	32(%r13),%r13
	leaq 	64(%r14),%r14
	addq 	$4,%rcx
	cmpq	%rax,%rcx
	jb 	fftoffloop
	/* end of block loop */
	leaq	(%rbx,%rax,2),%rbx
	cmpq	%r9,%rbx
	jb 	fftblockloop
	/* end of halfnn loop */
	shlq	$1,%rax
	leaq	(%rdx,%rax,8),%rdx
	cmpq	%r9,%rax
	jb ffthalfnnloop
















	/* cur_tt is at rdx */
	movq $0,%rax /* j */
	movq %rdi,%r10
	movq %rsi,%r11
fftfinalloop:
	vmovapd	(%r10),%ymm0 /* re */
	vmovapd	(%r11),%ymm1 /* im */
	vmovapd (%rdx),%ymm2 /* cos */
	vmovapd 32(%rdx),%ymm3 /* sin */
	vmulpd %ymm0,%ymm2,%ymm4
	vmulpd %ymm0,%ymm3,%ymm5
	vmulpd %ymm1,%ymm2,%ymm6
	vmulpd %ymm1,%ymm3,%ymm7
	vsubpd %ymm7,%ymm4,%ymm0
	vaddpd %ymm6,%ymm5,%ymm1
	vmovapd %ymm0,(%r10)
	vmovapd %ymm1,(%r11)
    	/* end of final loop */
    	leaq	32(%r10),%r10
    	leaq	32(%r11),%r11
    	leaq	64(%rdx),%rdx
	addq	$4,%rax
	cmpq	%r9,%rax
	jb fftfinalloop

	/* Restore registers */
fftend:
	vzeroall
	popq        %rbx
	popq        %r14
	popq        %r13
	popq        %r12
	popq        %r11
	popq        %r10
	retq


/* Constants for YMM */
.balign 32
size4negation0: .double +1.0, +1.0, +1.0, -1.0 /* ymm15 */
size4negation1: .double +1.0, -1.0, -1.0, +1.0 /* ymm14 */
size4negation2: .double +1.0, +1.0, -1.0, -1.0 /* ymm13 */
size4negation3: .double +1.0, -1.0, +1.0, -1.0 /* ymm12 */

	.cfi_endproc
.LFE0:
	.size	fft, .-fft
	.section	.text.unlikely
.LCOLDE0:
	.text
.LHOTE0:
	.ident	"GCC: (Ubuntu 5.2.1-22ubuntu2) 5.2.1 20151010"
	.section	.note.GNU-stack,"",@progbits

