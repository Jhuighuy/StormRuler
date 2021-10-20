// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< //
// Copyright (C) 2021 Oleg Butakov
// 
// Permission is hereby granted, free of charge, to any person 
// obtaining a copy of this software and associated documentation 
// files (the "Software"), to deal in the Software without 
// restriction, including without limitation the rights  to use, 
// copy, modify, merge, publish, distribute, sublicense, and/or
// sell copies of the Software, and to permit persons to whom the  
// Software is furnished to do so, subject to the following 
// conditions:
// 
// The above copyright notice and this permission notice shall be 
// included in all copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
// OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND 
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT 
// HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, 
// WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
// OTHER DEALINGS IN THE SOFTWARE.
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> //

#define _GNU_SOURCE

#include <unistd.h>
#include <sys/un.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netdb.h>
#include <netinet/in.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< //
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> //

/**
 * Read from the socket.
 * 
 * @param[in] sock Socket ID.
 * @param[out] data Data pointer.
 * @param[in] length Length of the data in bytes.
*/
void SR_SocketRead(int sock, char* data, size_t length) {
  while (length > 0) {
    ssize_t n = read(sock, data, length);
    if (n <= 0) {
      fputs("socket `read` failed.", stderr), abort();
    }
    data += n, length -= n;
  }
} // SR_SocketRead

/**
 * Write to the socket.
 *
 * @param[in] sock Socket ID.
 * @param[in] data Data pointer.
 * @param[in] length Length of the data in bytes.
*/
void SR_SocketWrite(int sock, const char* data, size_t length) {
  if (write(sock, data, length) < 0) {
    fputs("socket `write` failed.", stderr), abort();
  }
} // SR_SocketWrite

/**
 * Run a server.
 * 
 * @param[in] port Port of the server.
 * @param[in] callback Server callback. 
 */
void SR_RunServer(int port, void(*callback)(int fd)) {
	const int sock = socket(AF_INET, SOCK_STREAM, 0);
  if (sock < 0) {
    fputs("failed to create the socket.", stderr), abort();
  }

	struct sockaddr_in serv_addr = {};
	serv_addr.sin_family = AF_INET;
	serv_addr.sin_addr.s_addr = htonl(INADDR_ANY);
	serv_addr.sin_port = htons(port != 0 ? port : 2020);
	int error = bind(sock, (struct sockaddr*)&serv_addr, sizeof(serv_addr));
	error |= listen(sock, 1);
  if (error != 0) {
    fputs("failed to bind the socket.", stderr), abort();
  }

  const int connection = accept(sock, NULL, NULL);
  if (sock < 0) {
    fputs("failed to accept a request.", stderr), abort();
  }

  callback(connection);
  
  close(connection);
  close(sock);
} // SR_RunServer
