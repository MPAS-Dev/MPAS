#include <stdio.h>
#include <stdarg.h>

#define MAX_LINE_LEN 132

char printbuf[MAX_LINE_LEN+2];
char fbuffer[1024];
int nbuf = 0;

void fortprintf(FILE * fd, char * str, ...)
{
   int i, nl, sp, inquotes;
   int lastchar;
   va_list ap;

   /* Add formatted string to the buffer of fortran code to be written */
   va_start(ap, str);
   i = vsnprintf(fbuffer+nbuf, 1024-nbuf, str, ap);
   va_end(ap);

   /* Set the next free position in the fortran buffer */
   nbuf = nbuf + i;

   inquotes = 0;
   do {

      nl = sp = -1;

      /* Scan through the max line length - 1 (since we may have to add an & character) or the end of the buffer, whichever comes first */
      for(i=0; i<MAX_LINE_LEN-1 && i<nbuf; i++) {
         lastchar = (i == nbuf-1) ? 1 : 0;
         if (fbuffer[i] == '\'' && (fbuffer[i+1] != '\'' || lastchar)) inquotes = (inquotes + 1) % 2;  /* Whether we are inside a quoted string */
         if (fbuffer[i] == '\n') nl = i;                                                               /* The last occurrence of a newline */
         if (fbuffer[i] == ' ' && !lastchar && fbuffer[i+1] != '&') sp = i;                            /* The last occurrence of a space */
      }

      /* If the charater at column MAX_LINE_LEN happens to be a newline, though, we mark it */
      if (i == MAX_LINE_LEN-1 && fbuffer[i] == '\n') nl = i;

      /* If we haven't reached the column limit, don't consider breaking the line yet */
      if (nbuf <= MAX_LINE_LEN) sp = -1;

      /* If we have a newline */
      if (nl > 0) {
         snprintf(printbuf, nl+2, "%s", fbuffer);
         fprintf(fd, "%s", printbuf);
         nl++;

         /* Shift unprinted contents of fortran buffer to the beginning */
         for(i=0; nl<nbuf; i++, nl++)
            fbuffer[i] = fbuffer[nl];
         nbuf = i;
      }
      /* Else if we found a place to break the line */
      else if (sp > 0) {
         snprintf(printbuf, sp+2, "%s", fbuffer);
         i = sp+1;
         if (inquotes) printbuf[i++] = '\'';
         printbuf[i++] = '&';
         printbuf[i++] = '\n';
         printbuf[i++] = '\0';
         fprintf(fd, "%s", printbuf);
         sp++;
         i = 0;
         if (inquotes) {
            inquotes = (inquotes + 1) % 2;
            fbuffer[i++] = '/';
            fbuffer[i++] = '/';
            fbuffer[i++] = '\'';
         }

         /* Shift unprinted contents of fortran buffer to the beginning */
         for( ; sp<nbuf; i++, sp++)
            fbuffer[i] = fbuffer[sp];
         nbuf = i;
      }

   } while (nl > 0 || sp > 0);
}

void fortprint_flush(FILE * fd)
{
   snprintf(printbuf, nbuf+1, "%s", fbuffer);
   fprintf(fd, "%s", printbuf);
   nbuf = 0;
}
